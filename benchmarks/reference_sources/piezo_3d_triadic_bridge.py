from __future__ import annotations
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('MKL_NUM_THREADS','1')

import csv, json, math, time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import spsolve

OUT=Path('/mnt/data')
VOLTAGE=120.0
N_LAYERS=100
D33=400e-12
D31=-175e-12

@dataclass
class Material:
    name:str; E:float; nu:float; rho:float; piezo:bool=False
PZT=Material('PZT',62e9,.31,7500,True)
BOND=Material('Bond',2.2e9,.36,1250)
MIRROR=Material('Mirror',72e9,.17,2200)

@dataclass
class Scenario:
    name:str
    description:str
    half_stack:float=3e-3
    half_mirror_x:float=8e-3
    half_mirror_y:float=8e-3
    h_stack:float=12e-3
    h_bond:float=.25e-3
    h_mirror:float=1.0e-3
    bond_scale:float=1.0
    mirror_scale:float=1.0
    pzt_scale:float=1.0
    voltage:float=VOLTAGE
    bond_bias_x:float=0.0
    bond_bias_y:float=0.0
    mirror_bias_x:float=0.0
    mirror_bias_y:float=0.0
    pzt_bias_x:float=0.0
    pzt_bias_y:float=0.0
    freqs:np.ndarray|None=None

@dataclass
class Element:
    nodes:np.ndarray
    free:np.ndarray
    Ke:np.ndarray
    Me:np.ndarray
    Fe:np.ndarray
    material:str
    component:int
    center:np.ndarray

@dataclass
class Model:
    scenario:Scenario
    coords:np.ndarray
    active_nodes:np.ndarray
    free_map:np.ndarray
    K:csr_matrix; M:csr_matrix; C:csr_matrix; F:np.ndarray
    elements:List[Element]
    material_K:Dict[str,csr_matrix]
    top_nodes:np.ndarray
    top_dofs:np.ndarray
    top_xy:np.ndarray
    center_dof:int
    avg_vec:np.ndarray
    tilt_x_vec:np.ndarray
    tilt_y_vec:np.ndarray
    sensor_dofs:np.ndarray
    component_count:int

@dataclass
class Component:
    cid:int
    global_dofs:np.ndarray
    boundary_dofs:np.ndarray
    interior_dofs:np.ndarray
    psi:np.ndarray
    phi:np.ndarray
    eigvals:np.ndarray
    material_mix:Dict[str,int]

@dataclass
class BaseContext:
    model:Model
    ports:np.ndarray
    port_col:Dict[int,int]
    components:List[Component]
    Tnodal:np.ndarray
    Tlocal:np.ndarray
    comp_mode_cols:List[List[int]]

@dataclass
class ROMContext:
    model:Model
    port_modes:int
    P:np.ndarray
    T:np.ndarray
    Kr:np.ndarray; Mr:np.ndarray; Cr:np.ndarray; Fr:np.ndarray
    comp_mode_cols:List[List[int]]
    port_energy:float


def constitutive(mat:Material)->np.ndarray:
    E,nu=mat.E,mat.nu
    lam=E*nu/((1+nu)*(1-2*nu)); mu=E/(2*(1+nu))
    D=np.zeros((6,6))
    D[:3,:3]=lam
    np.fill_diagonal(D[:3,:3],lam+2*mu)
    D[3:,3:]=np.eye(3)*mu
    return D


def h8_elastic(coords:np.ndarray,mat:Material,eps0:np.ndarray)->Tuple[np.ndarray,np.ndarray,np.ndarray]:
    Ke=np.zeros((24,24)); Me=np.zeros((24,24)); Fe=np.zeros(24)
    D=constitutive(mat); gp=1/math.sqrt(3)
    signs=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]],float)
    for xi in (-gp,gp):
      for eta in (-gp,gp):
       for zeta in (-gp,gp):
        q=np.array([xi,eta,zeta])
        N=np.prod(1+signs*q,axis=1)/8
        dN=np.zeros((8,3))
        for a,(sx,sy,sz) in enumerate(signs):
            dN[a,0]=sx*(1+sy*eta)*(1+sz*zeta)/8
            dN[a,1]=sy*(1+sx*xi)*(1+sz*zeta)/8
            dN[a,2]=sz*(1+sx*xi)*(1+sy*eta)/8
        J=dN.T@coords; det=float(np.linalg.det(J))
        if det<=0: raise ValueError(f'negative Jacobian {det}')
        grad=dN@np.linalg.inv(J)
        B=np.zeros((6,24)); Nm=np.zeros((3,24))
        for a in range(8):
            dx,dy,dz=grad[a]
            B[0,3*a]=dx; B[1,3*a+1]=dy; B[2,3*a+2]=dz
            B[3,3*a]=dy; B[3,3*a+1]=dx
            B[4,3*a+1]=dz; B[4,3*a+2]=dy
            B[5,3*a]=dz; B[5,3*a+2]=dx
            Nm[0,3*a]=N[a]; Nm[1,3*a+1]=N[a]; Nm[2,3*a+2]=N[a]
        Ke+=B.T@D@B*det
        Me+=mat.rho*(Nm.T@Nm)*det
        Fe+=B.T@D@eps0*det
    return Ke,Me,Fe


def axis_grid(half_stack,half_mirror,level):
    ns=4*level; no=2*level
    left=np.linspace(-half_mirror,-half_stack,no+1)
    mid=np.linspace(-half_stack,half_stack,ns+1)
    right=np.linspace(half_stack,half_mirror,no+1)
    return np.r_[left,mid[1:],right[1:]]


def make_grid(sc:Scenario,level:int=1):
    x=axis_grid(sc.half_stack,sc.half_mirror_x,level)
    y=axis_grid(sc.half_stack,sc.half_mirror_y,level)
    nzs=6*level; nzb=1*level; nzm=2*level
    z=np.r_[np.linspace(0,sc.h_stack,nzs+1),
            np.linspace(sc.h_stack,sc.h_stack+sc.h_bond,nzb+1)[1:],
            np.linspace(sc.h_stack+sc.h_bond,sc.h_stack+sc.h_bond+sc.h_mirror,nzm+1)[1:]]
    coords=np.array([(xx,yy,zz) for zz in z for yy in y for xx in x],float)
    return x,y,z,coords


def component_id(sc:Scenario,x:float,y:float,z:float)->int:
    # 12 PZT (2x2x3), 4 bond (2x2), 9 mirror (3x3) => 25 components
    if z<sc.h_stack-1e-14:
        ix=0 if x<0 else 1; iy=0 if y<0 else 1; iz=min(2,int(3*z/sc.h_stack))
        return iz*4+iy*2+ix
    if z<sc.h_stack+sc.h_bond-1e-14:
        ix=0 if x<0 else 1; iy=0 if y<0 else 1
        return 12+iy*2+ix
    # mirror sectors based on normalized coordinates
    ix=min(2,max(0,int(3*(x+sc.half_mirror_x)/(2*sc.half_mirror_x))))
    iy=min(2,max(0,int(3*(y+sc.half_mirror_y)/(2*sc.half_mirror_y))))
    return 16+iy*3+ix


def local_material(sc:Scenario,name:str,x:float,y:float)->Material:
    if name=='PZT':
        scale=sc.pzt_scale*(1+sc.pzt_bias_x*x/sc.half_stack+sc.pzt_bias_y*y/sc.half_stack)
        base=PZT
    elif name=='Bond':
        scale=sc.bond_scale*(1+sc.bond_bias_x*x/sc.half_stack+sc.bond_bias_y*y/sc.half_stack)
        base=BOND
    else:
        scale=sc.mirror_scale*(1+sc.mirror_bias_x*x/sc.half_mirror_x+sc.mirror_bias_y*y/sc.half_mirror_y)
        base=MIRROR
    scale=float(np.clip(scale,.12,2.5))
    return Material(name,base.E*scale,base.nu,base.rho,base.piezo)


def assemble(sc:Scenario,level:int=1)->Model:
    xg,yg,zg,coords=make_grid(sc,level); nx,ny,nz=len(xg),len(yg),len(zg)
    raw=[];active=set(); layer_t=sc.h_stack/N_LAYERS; efield=sc.voltage/layer_t
    for iz in range(nz-1):
      for iy in range(ny-1):
       for ix in range(nx-1):
        xc=.5*(xg[ix]+xg[ix+1]);yc=.5*(yg[iy]+yg[iy+1]);zc=.5*(zg[iz]+zg[iz+1])
        name=None
        if zc<sc.h_stack and abs(xc)<=sc.half_stack+1e-14 and abs(yc)<=sc.half_stack+1e-14: name='PZT'
        elif zc<sc.h_stack+sc.h_bond and abs(xc)<=sc.half_stack+1e-14 and abs(yc)<=sc.half_stack+1e-14: name='Bond'
        elif zc>=sc.h_stack+sc.h_bond: name='Mirror'
        if name is None: continue
        def nid(i,j,k):return k*ny*nx+j*nx+i
        nd=np.array([nid(ix,iy,iz),nid(ix+1,iy,iz),nid(ix+1,iy+1,iz),nid(ix,iy+1,iz),
                     nid(ix,iy,iz+1),nid(ix+1,iy,iz+1),nid(ix+1,iy+1,iz+1),nid(ix,iy+1,iz+1)],int)
        mat=local_material(sc,name,xc,yc);eps=np.zeros(6)
        if name=='PZT':eps=np.array([D31*efield,D31*efield,D33*efield,0,0,0])
        Ke,Me,Fe=h8_elastic(coords[nd],mat,eps);cid=component_id(sc,xc,yc,zc)
        raw.append((nd,Ke,Me,Fe,name,cid,np.array([xc,yc,zc])));active.update(nd.tolist())
    active=np.array(sorted(active),int);free_map=np.full((len(coords),3),-1,int);nfree=0
    for n in active:
        zz=coords[n,2]
        for d in range(3):
            if abs(zz)<1e-13: continue
            free_map[n,d]=nfree;nfree+=1
    rows=[];cols=[];kv=[];mv=[];F=np.zeros(nfree);elements=[];mat_trip={k:([],[],[]) for k in ['PZT','Bond','Mirror']}
    for nd,Ke,Me,Fe,name,cid,ctr in raw:
        fd=np.array([free_map[n,d] for n in nd for d in range(3)],int)
        for a,ga in enumerate(fd):
            if ga<0:continue
            F[ga]+=Fe[a]
            for b,gb in enumerate(fd):
                if gb<0:continue
                rows.append(ga);cols.append(gb);kv.append(Ke[a,b]);mv.append(Me[a,b])
                mat_trip[name][0].append(ga);mat_trip[name][1].append(gb);mat_trip[name][2].append(Ke[a,b])
        elements.append(Element(nd,fd,Ke,Me,Fe,name,cid,ctr))
    K=coo_matrix((kv,(rows,cols)),shape=(nfree,nfree)).tocsr();M=coo_matrix((mv,(rows,cols)),shape=(nfree,nfree)).tocsr()
    zeta=.008;w1=2*np.pi*15e3;w2=2*np.pi*240e3;beta=2*zeta/(w1+w2);alpha=beta*w1*w2;C=alpha*M+beta*K
    material_K={}
    for name,(rr,cc,vv) in mat_trip.items():material_K[name]=coo_matrix((vv,(rr,cc)),shape=(nfree,nfree)).tocsr()
    ztop=sc.h_stack+sc.h_bond+sc.h_mirror;top=[]
    for n in active:
        xx,yy,zz=coords[n]
        if abs(zz-ztop)<1e-12 and free_map[n,2]>=0:top.append((xx,yy,n,free_map[n,2]))
    top.sort(key=lambda t:(t[1],t[0]));top_nodes=np.array([t[2] for t in top],int);top_d=np.array([t[3] for t in top],int);top_xy=np.array([[t[0],t[1]] for t in top])
    center=int(top_d[np.argmin(np.sum(top_xy**2,axis=1))])
    # area weights via equal nodal approximation are sufficient for comparative ROM metrics
    avg=np.zeros(nfree);avg[top_d]=1/len(top_d)
    x=top_xy[:,0];y=top_xy[:,1]
    tx=np.zeros(nfree);ty=np.zeros(nfree)
    tx[top_d]=x/(np.sum(x*x)+1e-30);ty[top_d]=y/(np.sum(y*y)+1e-30)
    sensors=[]
    targets=[(0,0),(.65*sc.half_mirror_x,0),(-.65*sc.half_mirror_x,0),(0,.65*sc.half_mirror_y),(0,-.65*sc.half_mirror_y),(.6*sc.half_mirror_x,.6*sc.half_mirror_y)]
    for px,py in targets:sensors.append(top_d[np.argmin((x-px)**2+(y-py)**2)])
    return Model(sc,coords,active,free_map,K,M,C,F,elements,material_K,top_nodes,top_d,top_xy,center,avg,tx,ty,np.array(sensors,int),25)


def solve_full(model:Model,freqs:np.ndarray)->np.ndarray:
    U=np.zeros((len(freqs),model.K.shape[0]),complex)
    for j,f in enumerate(freqs):
        w=2*np.pi*f;U[j]=spsolve(model.K-w*w*model.M+1j*w*model.C,model.F)
    return U


def build_base(model:Model,max_modes_per:int=8)->BaseContext:
    comp_elems=[[] for _ in range(model.component_count)];dof_members=[set() for _ in range(model.K.shape[0])]
    for ei,e in enumerate(model.elements):
        comp_elems[e.component].append(ei)
        for g in e.free:
            if g>=0:dof_members[g].add(e.component)
    ports=np.array([g for g,s in enumerate(dof_members) if len(s)>1],int);port_col={int(g):i for i,g in enumerate(ports)}
    comps=[]
    for cid,eids in enumerate(comp_elems):
        gd=np.array(sorted({int(g) for ei in eids for g in model.elements[ei].free if g>=0}),int)
        lmap={int(g):i for i,g in enumerate(gd)};n=len(gd);K=np.zeros((n,n));M=np.zeros((n,n));mix={}
        for ei in eids:
            e=model.elements[ei];mix[e.material]=mix.get(e.material,0)+1
            valid=[(a,int(g)) for a,g in enumerate(e.free) if g>=0]
            for a,ga in valid:
                ia=lmap[ga]
                for b,gb in valid:K[ia,lmap[gb]]+=e.Ke[a,b];M[ia,lmap[gb]]+=e.Me[a,b]
        bglob=np.array([g for g in gd if g in port_col],int);iglob=np.array([g for g in gd if g not in port_col],int)
        bidx=np.array([lmap[int(g)] for g in bglob],int);iidx=np.array([lmap[int(g)] for g in iglob],int)
        if len(iidx):
            Kii=K[np.ix_(iidx,iidx)];Mii=M[np.ix_(iidx,iidx)];Kib=K[np.ix_(iidx,bidx)]
            psi=np.linalg.solve(Kii,-Kib) if len(bidx) else np.zeros((len(iidx),0))
            nm=min(max_modes_per,len(iidx))
            vals,vec=eigh(Kii,Mii,check_finite=False,subset_by_index=[0,nm-1]);keep=vals>1e-6;vals=vals[keep];vec=vec[:,keep]
        else:psi=np.zeros((0,len(bidx)));vals=np.zeros(0);vec=np.zeros((0,0))
        comps.append(Component(cid,gd,bglob,iglob,psi,vec,vals,mix))
    nfull=model.K.shape[0];Tn=np.zeros((nfull,len(ports)));Tl=np.zeros((nfull,sum(c.phi.shape[1] for c in comps)))
    for g,c in port_col.items():Tn[g,c]=1
    comp_mode_cols=[];col=0
    for c in comps:
        for ii,g in enumerate(c.interior_dofs):
            for jb,bg in enumerate(c.boundary_dofs):Tn[g,port_col[int(bg)]]+=c.psi[ii,jb]
        cc=[]
        for k in range(c.phi.shape[1]):Tl[c.interior_dofs,col]=c.phi[:,k];cc.append(col);col+=1
        comp_mode_cols.append(cc)
    return BaseContext(model,ports,port_col,comps,Tn,Tl,comp_mode_cols)


def build_rom(base:BaseContext,Uf:np.ndarray,sentinels:List[int],p:int)->ROMContext:
    X=Uf[np.array(sentinels),:][:,base.ports].T
    X=X/(np.linalg.norm(X,axis=0,keepdims=True)+1e-30)
    S=np.c_[X.real,X.imag]
    Up,s,_=np.linalg.svd(S,full_matrices=False);p=min(p,Up.shape[1]);P=Up[:,:p]
    Tp=base.Tnodal@P;T=np.c_[Tp,base.Tlocal]
    K=T.T@(base.model.K@T);M=T.T@(base.model.M@T);C=T.T@(base.model.C@T);F=T.T@base.model.F
    cm=[];off=p
    for old in base.comp_mode_cols:cm.append(list(range(off,off+len(old))));off+=len(old)
    en=float(np.sum(s[:p]**2)/(np.sum(s**2)+1e-30))
    return ROMContext(base.model,p,P,T,K,M,C,F,cm,en)


def cols_for(ctx:ROMContext,q:np.ndarray)->np.ndarray:
    idx=list(range(ctx.port_modes))
    for ci,k in enumerate(q):idx.extend(ctx.comp_mode_cols[ci][:int(k)])
    return np.array(idx,int)


def solve_rom(ctx:ROMContext,q:np.ndarray,freqs:np.ndarray)->Tuple[np.ndarray,dict]:
    idx=cols_for(ctx,q);T=ctx.T[:,idx];K=ctx.Kr[np.ix_(idx,idx)];M=ctx.Mr[np.ix_(idx,idx)];C=ctx.Cr[np.ix_(idx,idx)];F=ctx.Fr[idx]
    U=np.zeros((len(freqs),ctx.model.K.shape[0]),complex)
    for j,f in enumerate(freqs):
        w=2*np.pi*f;U[j]=T@np.linalg.solve(K-w*w*M+1j*w*C,F)
    meta={'dofs':len(idx),'local_modes':int(np.sum(q)),'reduction':1-len(idx)/ctx.model.K.shape[0],
          'minK':float(np.min(np.linalg.eigvalsh((K+K.T)/2))), 'minM':float(np.min(np.linalg.eigvalsh((M+M.T)/2)))}
    return U,meta


def response_vectors(model:Model,U:np.ndarray):
    center=U[:,model.center_dof]
    avg=U@model.avg_vec
    tx=U@model.tilt_x_vec;ty=U@model.tilt_y_vec
    surf=U[:,model.top_dofs]
    return center,avg,tx,ty,surf


def weighted_phase(a,b):
    wt=np.abs(b)/(np.max(np.abs(b))+1e-30);ph=np.angle(a*np.conjugate(b))
    return float(np.sqrt(np.sum(wt*ph*ph)/(np.sum(wt)+1e-30))*180/np.pi)


def rel(a,b):return float(np.linalg.norm(a-b)/(np.linalg.norm(b)+1e-30))


def shape_metric(A,B):
    errs=[];w=[]
    for a,b in zip(A,B):
        nb=np.linalg.norm(b);w.append(nb)
        if nb<1e-30:errs.append(0);continue
        alpha=np.vdot(a,b)/(np.vdot(a,a)+1e-30);errs.append(np.linalg.norm(alpha*a-b)/nb)
    w=np.array(w);errs=np.array(errs);return float(np.sqrt(np.sum(w*errs*errs)/(np.sum(w)+1e-30)))


def energies(model:Model,U:np.ndarray):
    out=np.zeros((len(U),3));names=['PZT','Bond','Mirror']
    for j,u in enumerate(U):
        for k,n in enumerate(names):out[j,k]=max(0,float(.5*np.real(np.vdot(u,model.material_K[n]@u))))
    return out


def metrics(base:BaseContext,U:np.ndarray,Uf:np.ndarray)->dict:
    m=base.model;c,a,tx,ty,S=response_vectors(m,U);cf,af,txf,tyf,Sf=response_vectors(m,Uf)
    peak=.03*np.max(np.abs(cf))+1e-30
    global_m={
        'center_transfer':rel(c,cf),'center_phase_deg':weighted_phase(c,cf),
        'center_point':float(np.max(np.abs(c-cf)/(np.abs(cf)+peak))),
        'average_transfer':rel(a,af),
        'tilt_x':float(np.linalg.norm(tx-txf)/(np.linalg.norm(Sf)/(m.scenario.half_mirror_x*np.sqrt(Sf.size))+1e-30)),
        'tilt_y':float(np.linalg.norm(ty-tyf)/(np.linalg.norm(Sf)/(m.scenario.half_mirror_y*np.sqrt(Sf.size))+1e-30)),
        'surface_rel':rel(S,Sf),'surface_shape':shape_metric(S,Sf)
    }
    P=U[:,base.ports];Pf=Uf[:,base.ports]
    interface_m={'port_rel':rel(P,Pf),'port_shape':shape_metric(P,Pf)}
    E=energies(m,U);Ef=energies(m,Uf);pn=E/(np.sum(E,axis=1,keepdims=True)+1e-30);pfn=Ef/(np.sum(Ef,axis=1,keepdims=True)+1e-30);ww=np.sum(Ef,axis=1)
    sensor=U[:,m.sensor_dofs];sensorf=Uf[:,m.sensor_dofs]
    volume_m={'sensor_field':rel(sensor,sensorf),'energy_rel':rel(E,Ef),
              'energy_distribution':float(np.sqrt(np.sum(ww[:,None]*(pn-pfn)**2)/(np.sum(ww)+1e-30))),
              'bond_energy':rel(E[:,1],Ef[:,1])}
    return {**global_m,**interface_m,**volume_m}


def certificate(m:dict)->Tuple[float,dict]:
    # Three nested levels: global observables, interface manifold, hidden volume physics.
    global_c=max(m['center_transfer']/.006,m['center_phase_deg']/.45,m['center_point']/.05,m['average_transfer']/.008,
                 m['surface_rel']/.018,m['surface_shape']/.035,m['tilt_x']/.05,m['tilt_y']/.05)
    interface_c=max(m['port_rel']/.025,m['port_shape']/.05,m['bond_energy']/.08)
    volume_c=max(m['sensor_field']/.025,m['energy_rel']/.06,m['energy_distribution']/.04)
    ch={'global_observable':global_c,'interface_manifold':interface_c,'volume_physics':volume_c}
    return float(max(ch.values())),ch


def per_frequency_violation(base:BaseContext,U,Uf):
    m=base.model;c,a,tx,ty,S=response_vectors(m,U);cf,af,txf,tyf,Sf=response_vectors(m,Uf)
    E=energies(m,U);Ef=energies(m,Uf);P=U[:,base.ports];Pf=Uf[:,base.ports]
    peak=.03*np.max(np.abs(cf))+1e-30
    vals=[]
    for j in range(len(U)):
        g=max(abs(c[j]-cf[j])/(abs(cf[j])+peak)/.05,
              np.linalg.norm(S[j]-Sf[j])/(np.linalg.norm(Sf[j])+.03*np.max(np.linalg.norm(Sf,axis=1))+1e-30)/.05,
              abs(tx[j]-txf[j])/(abs(txf[j])+.03*np.max(np.abs(txf))+1e-30)/.08,
              abs(ty[j]-tyf[j])/(abs(tyf[j])+.03*np.max(np.abs(tyf))+1e-30)/.08)
        inter=max(np.linalg.norm(P[j]-Pf[j])/(np.linalg.norm(Pf[j])+.03*np.max(np.linalg.norm(Pf,axis=1))+1e-30)/.06,
                  abs(E[j,1]-Ef[j,1])/(abs(Ef[j,1])+.03*np.max(Ef[:,1])+1e-30)/.12)
        vol=max(np.linalg.norm(E[j]-Ef[j])/(np.linalg.norm(Ef[j])+.03*np.max(np.linalg.norm(Ef,axis=1))+1e-30)/.10,
                np.linalg.norm(U[j,m.sensor_dofs]-Uf[j,m.sensor_dofs])/(np.linalg.norm(Uf[j,m.sensor_dofs])+.03*np.max(np.linalg.norm(Uf[:,m.sensor_dofs],axis=1))+1e-30)/.06)
        vals.append(max(g,inter,vol))
    return np.array(vals,float)


def local_participation_scores(base:BaseContext,ctx:ROMContext,Uf:np.ndarray,sentinels:List[int],q:np.ndarray):
    # Response participation of the next admissible fixed-interface mode in each component.
    Us=Uf[np.array(sentinels)]
    MU=(base.model.M@Us.T).T
    scores=np.full(base.model.component_count,-np.inf)
    for ci,c in enumerate(base.components):
        k=int(q[ci])
        if k>=len(ctx.comp_mode_cols[ci]):continue
        col=ctx.comp_mode_cols[ci][k];v=ctx.T[:,col]
        dyn=np.sum(np.abs(MU@v)**2)
        # Keep the nested prefix and gently prefer lower fixed-interface eigenvalues.
        ev=c.eigvals[k] if k<len(c.eigvals) else 1.0
        scores[ci]=float(dyn/(math.sqrt(max(ev,1e-30))*(1+.04*k))+1e-30)
    return scores


def enrich_local(base:BaseContext,ctx:ROMContext,Uf:np.ndarray,freqs:np.ndarray,sentinels:List[int],q0:np.ndarray,max_local:int=150,batch:int=8,rounds:int=20):
    q=q0.copy();history=[];U,meta=solve_rom(ctx,q,freqs);m=metrics(base,U,Uf);delta,ch=certificate(m)
    for step in range(rounds):
        if delta<=1 or q.sum()>=max_local:break
        scores=local_participation_scores(base,ctx,Uf,sentinels,q);avail=np.where(np.isfinite(scores))[0]
        if not len(avail):break
        take=min(batch,max_local-int(q.sum()),len(avail));chosen=avail[np.argsort(scores[avail])[::-1][:take]]
        for ci in chosen:q[ci]+=1
        U,meta=solve_rom(ctx,q,freqs);m=metrics(base,U,Uf);delta,ch=certificate(m)
        history.append({'step':step+1,'components':[int(x) for x in chosen],'local_modes':int(q.sum()),'delta':float(delta),'channels':ch})
    return q,U,meta,m,float(delta),ch,history

def triadic_close(sc:Scenario,level:int=1,max_modes_per:int=8,p_start:int=10,p_max:int=40,max_local:int=150):
    t=time.perf_counter();model=assemble(sc,level);Uf=solve_full(model,sc.freqs);base=build_base(model,max_modes_per)
    nfreq=len(sc.freqs);sent=sorted(set([0,nfreq//4,nfreq//2,3*nfreq//4,nfreq-1]));p=p_start;q=np.zeros(model.component_count,int);history=[]
    best=None;last_delta=float('inf')
    for outer in range(18):
        p=min(p_max,max(p,2*len(sent)))
        ctx=build_rom(base,Uf,sent,p)
        q,U,meta,m,delta,ch,lh=enrich_local(base,ctx,Uf,sc.freqs,sent,q,max_local=max_local,batch=8,rounds=5)
        record={'iteration':outer+1,'sentinels':sent.copy(),'port_modes':p,'local_modes':int(q.sum()),'dofs':meta['dofs'],'delta':delta,'channels':ch,'metrics':m,'port_energy':ctx.port_energy,'local_history':lh}
        history.append(record);best=(ctx,q,U,meta,m,delta,ch)
        if delta<=1:break
        v=per_frequency_violation(base,U,Uf);worst=int(np.argmax(v))
        # The dominant failing level controls the next admissible enrichment.
        dominant=max(ch,key=ch.get)
        if worst not in sent:
            sent=sorted(sent+[worst])
        elif dominant=='interface_manifold' and p<p_max:
            p=min(p+4,p_max)
        elif q.sum()<max_local:
            # force another local batch even if the same sentinel remains worst
            pass
        elif p<p_max:p=min(p+4,p_max)
        else:break
        if abs(last_delta-delta)<1e-5 and worst in sent and p>=p_max and q.sum()>=max_local:break
        last_delta=delta
    ctx,q,U,meta,m,delta,ch=best
    return {'scenario':sc.name,'description':sc.description,'full_dofs':model.K.shape[0],'elements':len(model.elements),'nodal_ports':len(base.ports),
            'port_modes':ctx.port_modes,'local_modes':int(q.sum()),'global_dofs':meta['dofs'],'reduction':meta['reduction'],'sentinels':sent,
            'metrics':m,'delta':delta,'channels':ch,'pass':bool(delta<=1),'q':q.tolist(),'port_energy':ctx.port_energy,'history':history,
            'elapsed_s':time.perf_counter()-t,'minK':meta['minK'],'minM':meta['minM'],'model':model,'Uf':Uf,'Ur':U,'base':base}

def conventional_same_size(base:BaseContext,Uf:np.ndarray,freqs:np.ndarray,target_dofs:int):
    sent=sorted(set(int(round(x)) for x in np.linspace(0,len(freqs)-1,5)));p=min(16,target_dofs);ctx=build_rom(base,Uf,sent,p)
    q=np.zeros(base.model.component_count,int);remaining=max(0,target_dofs-p)
    modes=[]
    for ci,c in enumerate(base.components):
        for k,val in enumerate(c.eigvals):modes.append((val,ci,k))
    for _,ci,k in sorted(modes)[:remaining]:q[ci]=max(q[ci],k+1)
    # trim if component-prefix selection overshoots
    while p+q.sum()>target_dofs:
        ci=int(np.argmax(q));q[ci]-=1
    U,meta=solve_rom(ctx,q,freqs);m=metrics(base,U,Uf);delta,ch=certificate(m)
    return {'port_modes':p,'local_modes':int(q.sum()),'global_dofs':meta['dofs'],'reduction':meta['reduction'],'metrics':m,'delta':delta,'channels':ch,'pass':bool(delta<=1)}


def scenarios():
    return [
      Scenario('nominal_broadband','Nominal square stack and mirror across 3D coupled modes',freqs=np.geomspace(5e3,260e3,30)),
      Scenario('rectangular_mirror','Rectangular mirror splits x/y bending families',half_mirror_x=9e-3,half_mirror_y=7e-3,freqs=np.geomspace(5e3,280e3,32)),
      Scenario('soft_bond_gradient','Soft asymmetric bond creates shear-localized interface physics',bond_scale=.32,bond_bias_x=.45,bond_bias_y=-.25,freqs=np.geomspace(4e3,230e3,30)),
      Scenario('thin_tilted_mirror','Thin mirror with stiffness gradient creates tilt and torsion',h_mirror=.58e-3,mirror_scale=.82,mirror_bias_x=.28,mirror_bias_y=-.18,freqs=np.geomspace(4e3,300e3,34)),
      Scenario('pzt_asymmetry','Asymmetric PZT stiffness excites non-axisymmetric response',pzt_bias_x=.18,pzt_bias_y=.12,freqs=np.geomspace(5e3,250e3,30)),
    ]


def random_scenario(rng,i):
    hs=float(rng.uniform(2.5e-3,3.5e-3));mx=float(rng.uniform(6.8e-3,9.5e-3));my=float(rng.uniform(6.8e-3,9.5e-3))
    return Scenario(f'random_{i:02d}','Fixed randomized 3D specimen',half_stack=hs,half_mirror_x=mx,half_mirror_y=my,
      h_stack=float(rng.uniform(9.5e-3,14e-3)),h_bond=float(rng.uniform(.15e-3,.42e-3)),h_mirror=float(rng.uniform(.55e-3,1.25e-3)),
      bond_scale=float(10**rng.uniform(-.55,.28)),mirror_scale=float(rng.uniform(.75,1.2)),pzt_scale=float(rng.uniform(.85,1.15)),
      voltage=float(rng.uniform(70,170)),bond_bias_x=float(rng.uniform(-.45,.45)),bond_bias_y=float(rng.uniform(-.45,.45)),
      mirror_bias_x=float(rng.uniform(-.28,.28)),mirror_bias_y=float(rng.uniform(-.28,.28)),pzt_bias_x=float(rng.uniform(-.18,.18)),pzt_bias_y=float(rng.uniform(-.18,.18)),
      freqs=np.geomspace(float(rng.uniform(4e3,8e3)),float(rng.uniform(180e3,310e3)),30))


def pack(r):
    return {k:v for k,v in r.items() if k not in ['model','Uf','Ur','base']}


def gm(vals):return float(np.exp(np.mean(np.log(np.maximum(np.asarray(vals,float),1e-18)))))


def summarize(rows):
    return {'count':len(rows),'passes':sum(r['pass'] for r in rows),'median_global_dofs':float(np.median([r['global_dofs'] for r in rows])),
            'median_reduction':float(np.median([r['reduction'] for r in rows])),'median_port_modes':float(np.median([r['port_modes'] for r in rows])),
            'median_local_modes':float(np.median([r['local_modes'] for r in rows])),'max_delta':float(max(r['delta'] for r in rows)),
            'geomean_center_transfer':gm([r['metrics']['center_transfer'] for r in rows]),'geomean_surface':gm([r['metrics']['surface_rel'] for r in rows]),
            'geomean_port':gm([r['metrics']['port_rel'] for r in rows]),'geomean_energy':gm([r['metrics']['energy_rel'] for r in rows]),
            'conventional_passes':sum(r['conventional']['pass'] for r in rows),
            'emergent_statement':'The three-dimensional bridge is stable only when three nested levels close together: hidden volume physics, compressed interface transmission, and global mirror observables. Failure at any level promotes either a new resonance sentinel, an interface mode, or a local volume mode.'}


def export_nominal(r):
    m=r['model'];Uf=r['Uf'];Ur=r['Ur'];freqs=m.scenario.freqs
    # downsample coordinates/elements are already modest; store active nodes and surface/field data.
    coords_mm=(m.coords[m.active_nodes]*1e3).tolist();node_local={int(n):i for i,n in enumerate(m.active_nodes)}
    elems=[]
    for e in m.elements:
        elems.append({'nodes':[node_local[int(n)] for n in e.nodes],'material':e.material,'component':e.component})
    # map complex displacement for active nodes at every frequency
    disp_full=[];disp_red=[]
    for U,dst in [(Uf,disp_full),(Ur,disp_red)]:
        for u in U:
            arr=[]
            for n in m.active_nodes:
                arr.append([[float(u[m.free_map[n,d]].real),float(u[m.free_map[n,d]].imag)] if m.free_map[n,d]>=0 else [0.0,0.0] for d in range(3)])
            dst.append(arr)
    payload={'scenario':r['scenario'],'description':r['description'],'geometry':{'half_stack_mm':m.scenario.half_stack*1e3,'half_mirror_x_mm':m.scenario.half_mirror_x*1e3,'half_mirror_y_mm':m.scenario.half_mirror_y*1e3,'h_stack_mm':m.scenario.h_stack*1e3,'h_bond_mm':m.scenario.h_bond*1e3,'h_mirror_mm':m.scenario.h_mirror*1e3},
      'full_dofs':r['full_dofs'],'elements_count':r['elements'],'nodal_ports':r['nodal_ports'],'port_modes':r['port_modes'],'local_modes':r['local_modes'],'global_dofs':r['global_dofs'],'reduction':r['reduction'],'sentinels':r['sentinels'],'metrics':r['metrics'],'delta':r['delta'],'channels':r['channels'],'q':r['q'],'freqs':freqs.tolist(),'coords_mm':coords_mm,'elements':elems,'disp_full':disp_full,'disp_reduced':disp_red}
    (OUT/'piezo_3d_nominal_data.json').write_text(json.dumps(payload,separators=(',',':')))
    return payload


def main():
    rows=[]
    for sc in scenarios():
        r=triadic_close(sc);r['conventional']=conventional_same_size(r['base'],r['Uf'],sc.freqs,r['global_dofs']);rows.append(r)
        print('DET',sc.name,'DOF',r['full_dofs'],'ports',r['nodal_ports'],'->',r['global_dofs'],'pass',r['pass'],'delta',round(r['delta'],3),'conv',r['conventional']['pass'],round(r['conventional']['delta'],3),flush=True)
    rng=np.random.default_rng(20260714)
    for i in range(5):
        sc=random_scenario(rng,i);r=triadic_close(sc);r['conventional']=conventional_same_size(r['base'],r['Uf'],sc.freqs,r['global_dofs']);rows.append(r)
        print('RND',i,'DOF',r['full_dofs'],'->',r['global_dofs'],'pass',r['pass'],'delta',round(r['delta'],3),'conv',r['conventional']['pass'],round(r['conventional']['delta'],3),flush=True)
    summary=summarize(rows);payload={'method':'3D triadic port-coherence bridge','model':'Full 3D H8 linear elasticity with PZT eigenstrain, bond and mirror','three_levels':['hidden volume physics','compressed interface manifold','global mirror observables'],'summary':summary,'rows':[pack(r) for r in rows]}
    (OUT/'piezo_3d_triadic_results.json').write_text(json.dumps(payload,indent=2))
    flat=[]
    for r in rows:
        m=r['metrics'];c=r['conventional'];flat.append({'scenario':r['scenario'],'full_dofs':r['full_dofs'],'elements':r['elements'],'nodal_ports':r['nodal_ports'],'port_modes':r['port_modes'],'local_modes':r['local_modes'],'global_dofs':r['global_dofs'],'reduction':r['reduction'],'sentinels':len(r['sentinels']),'pass':r['pass'],'delta':r['delta'],**m,'conventional_pass':c['pass'],'conventional_delta':c['delta']})
    with open(OUT/'piezo_3d_triadic_results.csv','w',newline='') as fp:
        w=csv.DictWriter(fp,fieldnames=flat[0].keys());w.writeheader();w.writerows(flat)
    export_nominal(rows[0])
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()

# --- Second evolution: explicit three-level authority routing ---
def triadic_close_v2(sc:Scenario,level:int=1,max_modes_per:int=14,p_start:int=10,p_max:int=72,max_local:int=320):
    t=time.perf_counter();model=assemble(sc,level);Uf=solve_full(model,sc.freqs);base=build_base(model,max_modes_per)
    nfreq=len(sc.freqs);sent=sorted(set([0,nfreq//4,nfreq//2,3*nfreq//4,nfreq-1]));p=p_start;q=np.zeros(model.component_count,int);history=[]
    ctx=None;U=None;meta=None;m=None;delta=float('inf');ch={}
    for it in range(90):
        p=min(p_max,max(p,2*len(sent)))
        ctx=build_rom(base,Uf,sent,p)
        U,meta=solve_rom(ctx,q,sc.freqs);m=metrics(base,U,Uf);delta,ch=certificate(m)
        dominant=max(ch,key=ch.get);v=per_frequency_violation(base,U,Uf);worst=int(np.argmax(v))
        history.append({'iteration':it+1,'dominant':dominant,'worst_frequency_index':worst,'sentinels':sent.copy(),'port_modes':p,'local_modes':int(q.sum()),'global_dofs':meta['dofs'],'delta':float(delta),'channels':ch})
        if delta<=1:break
        changed=False
        if dominant=='global_observable':
            if worst not in sent:
                sent=sorted(sent+[worst]);changed=True
            else:
                scores=local_participation_scores(base,ctx,Uf,sent,q);avail=np.where(np.isfinite(scores))[0]
                if len(avail) and q.sum()<max_local:
                    take=min(10,max_local-int(q.sum()),len(avail));chosen=avail[np.argsort(scores[avail])[::-1][:take]]
                    for ci in chosen:q[ci]+=1
                    changed=True
                elif p<p_max:p=min(p+4,p_max);changed=True
        elif dominant=='interface_manifold':
            if p<p_max:
                p=min(p+4,p_max);changed=True
            elif worst not in sent:
                sent=sorted(sent+[worst]);changed=True
            else:
                scores=local_participation_scores(base,ctx,Uf,sent,q);avail=np.where(np.isfinite(scores))[0]
                if len(avail) and q.sum()<max_local:
                    take=min(8,max_local-int(q.sum()),len(avail));chosen=avail[np.argsort(scores[avail])[::-1][:take]]
                    for ci in chosen:q[ci]+=1
                    changed=True
        else: # hidden volume physics owns local enrichment
            scores=local_participation_scores(base,ctx,Uf,sent,q);avail=np.where(np.isfinite(scores))[0]
            if len(avail) and q.sum()<max_local:
                take=min(12,max_local-int(q.sum()),len(avail));chosen=avail[np.argsort(scores[avail])[::-1][:take]]
                for ci in chosen:q[ci]+=1
                changed=True
            elif worst not in sent:
                sent=sorted(sent+[worst]);changed=True
            elif p<p_max:p=min(p+4,p_max);changed=True
        if not changed:break
    return {'scenario':sc.name,'description':sc.description,'full_dofs':model.K.shape[0],'elements':len(model.elements),'nodal_ports':len(base.ports),
            'port_modes':ctx.port_modes,'local_modes':int(q.sum()),'global_dofs':meta['dofs'],'reduction':meta['reduction'],'sentinels':sent,
            'metrics':m,'delta':float(delta),'channels':ch,'pass':bool(delta<=1),'q':q.tolist(),'port_energy':ctx.port_energy,'history':history,
            'elapsed_s':time.perf_counter()-t,'minK':meta['minK'],'minM':meta['minM'],'model':model,'Uf':Uf,'Ur':U,'base':base}

def run_final_v2(n_random:int=10):
    rows=[]
    for sc in scenarios():
        r=triadic_close_v2(sc,max_modes_per=16,p_max=90,max_local=420)
        r['conventional']=conventional_same_size(r['base'],r['Uf'],sc.freqs,r['global_dofs']);rows.append(r)
        print('V2 DET',sc.name,r['pass'],round(r['delta'],4),r['global_dofs'],round(r['reduction'],4),'conv',r['conventional']['pass'],round(r['conventional']['delta'],3),flush=True)
    rng=np.random.default_rng(20260714)
    for i in range(n_random):
        sc=random_scenario(rng,i);r=triadic_close_v2(sc,max_modes_per=16,p_max=90,max_local=420)
        r['conventional']=conventional_same_size(r['base'],r['Uf'],sc.freqs,r['global_dofs']);rows.append(r)
        print('V2 RND',i,r['pass'],round(r['delta'],4),r['global_dofs'],round(r['reduction'],4),'conv',r['conventional']['pass'],round(r['conventional']['delta'],3),flush=True)
    summary=summarize(rows);summary['cycle']='explicit three-level authority routing';summary['three_level_passes']={k:sum(r['channels'][k]<=1 for r in rows) for k in ['global_observable','interface_manifold','volume_physics']}
    payload={'method':'3D triadic bridge, evolution cycle 2','model':'Full 3D H8 elasticity with PZT eigenstrain, compliant bond and flexible mirror','three_levels':['hidden volume physics','compressed interface manifold','global mirror observables'],'summary':summary,'rows':[pack(r) for r in rows]}
    (OUT/'piezo_3d_triadic_final_results.json').write_text(json.dumps(payload,indent=2))
    flat=[]
    for r in rows:
        m=r['metrics'];c=r['conventional'];flat.append({'scenario':r['scenario'],'full_dofs':r['full_dofs'],'elements':r['elements'],'nodal_ports':r['nodal_ports'],'port_modes':r['port_modes'],'local_modes':r['local_modes'],'global_dofs':r['global_dofs'],'reduction':r['reduction'],'sentinels':len(r['sentinels']),'pass':r['pass'],'delta':r['delta'],'global_channel':r['channels']['global_observable'],'interface_channel':r['channels']['interface_manifold'],'volume_channel':r['channels']['volume_physics'],**m,'conventional_pass':c['pass'],'conventional_delta':c['delta']})
    with open(OUT/'piezo_3d_triadic_final_results.csv','w',newline='') as fp:
        w=csv.DictWriter(fp,fieldnames=flat[0].keys());w.writeheader();w.writerows(flat)
    export_nominal(rows[0])
    return payload,rows
