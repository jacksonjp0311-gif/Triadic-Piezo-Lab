import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
from dataclasses import dataclass

L = 0.020
A0 = 16e-6
RHO0 = 7600.0
E0 = 62e9
D33 = 400e-12
VOLTAGE = 120.0
N_LAYERS = 120
MIRROR_MASS = 0.0018
MIRROR_K = 2.2e7
ZETA = 0.012


def Efun(x):
    x = np.asarray(x)
    bond = 0.58*np.exp(-((x-0.0142)/0.00062)**2)
    interface = 0.28*np.exp(-((x-0.0180)/0.00038)**2)
    ripple = 0.055*np.sin(2*np.pi*x/0.0025)**2
    return E0*np.clip(1.0 - bond - interface - ripple, 0.24, None)

def Afun(x):
    x = np.asarray(x)
    neck = 0.42*np.exp(-((x-0.0185)/0.00058)**2)
    return A0*(1.0 - neck)

def rhofun(x):
    return RHO0*np.ones_like(np.asarray(x))

@dataclass
class Model:
    x: np.ndarray
    K: csr_matrix
    M: csr_matrix
    C: csr_matrix
    F: np.ndarray


def assemble(x: np.ndarray) -> Model:
    n = len(x)
    K = lil_matrix((n,n), dtype=float)
    M = lil_matrix((n,n), dtype=float)
    F = np.zeros(n, float)
    eps = D33 * VOLTAGE * N_LAYERS / L
    for e in range(n-1):
        a,b = x[e],x[e+1]
        h = b-a
        xm = 0.5*(a+b)
        E = float(Efun(xm)); A=float(Afun(xm)); rho=float(rhofun(xm))
        c=E*A/h
        ke=((c,-c),(-c,c))
        mc=rho*A*h/6
        me=((2*mc,mc),(mc,2*mc))
        fe=E*A*eps
        i=e; j=e+1
        K[i,i]+=ke[0][0]; K[i,j]+=ke[0][1]; K[j,i]+=ke[1][0]; K[j,j]+=ke[1][1]
        M[i,i]+=me[0][0]; M[i,j]+=me[0][1]; M[j,i]+=me[1][0]; M[j,j]+=me[1][1]
        F[i]-=fe; F[j]+=fe
    M[-1,-1] += MIRROR_MASS
    K[-1,-1] += MIRROR_K
    K=K.tocsr(); M=M.tocsr()
    w1=2*np.pi*28e3; w2=2*np.pi*210e3
    beta = 2*ZETA/(w1+w2)
    alpha = beta*w1*w2
    C=(alpha*M + beta*K).tocsr()
    return Model(x,K,M,C,F)


def solve_harmonic(model: Model, freqs: np.ndarray) -> np.ndarray:
    K=model.K[1:,1:]; M=model.M[1:,1:]; C=model.C[1:,1:]; F=model.F[1:]
    U=np.zeros((len(freqs),len(model.x)),complex)
    for j,f in enumerate(freqs):
        w=2*np.pi*f
        Z=(K-w*w*M+1j*w*C).tocsc()
        U[j,1:]=spsolve(Z,F)
    return U


def recovered_stress_indicator(model: Model, U: np.ndarray, freqs: np.ndarray, weights=None) -> np.ndarray:
    x=model.x; ne=len(x)-1; nf=U.shape[0]
    if weights is None: weights=np.ones(nf)
    eta=np.zeros(ne,float)
    h=np.diff(x); mid=.5*(x[:-1]+x[1:]); E=Efun(mid); A=Afun(mid); rho=rhofun(mid)
    for j in range(nf):
        u=U[j]
        sig=E*np.diff(u)/h
        sn=np.zeros(len(x),complex)
        sn[0]=sig[0]; sn[-1]=sig[-1]
        for i in range(1,len(x)-1):
            wl=h[i-1]; wr=h[i]
            sn[i]=(wr*sig[i-1]+wl*sig[i])/(wl+wr)
        diffL=sn[:-1]-sig; diffR=sn[1:]-sig
        zz=h*A/np.maximum(E,1.0)*0.5*(np.abs(diffL)**2+np.abs(diffR)**2)
        um=.5*(u[:-1]+u[1:])
        w=2*np.pi*freqs[j]
        dyn=np.abs(rho*A*w*w*um)
        dr=(h**3/np.maximum(E*A,1e-30))*dyn**2/12
        eta += float(weights[j])*(zz+0.12*dr)
    return np.sqrt(np.maximum(eta,0))


def static_indicator(model: Model, u: np.ndarray) -> np.ndarray:
    return recovered_stress_indicator(model,u[None,:],np.array([0.0]),np.ones(1))


def interpolate_complex(xc, uc, xr):
    return np.interp(xr,xc,uc.real)+1j*np.interp(xr,xc,uc.imag)


def errors(model: Model, U: np.ndarray, xref: np.ndarray, Uref: np.ndarray, freqs: np.ndarray):
    tip=U[:,-1]; tipref=Uref[:,-1]
    qoi=np.linalg.norm(tip-tipref)/np.linalg.norm(tipref)
    num=0.; den=0.
    for j in range(len(freqs)):
        ui=interpolate_complex(model.x,U[j],xref)
        num += np.trapezoid(np.abs(ui-Uref[j])**2,xref)
        den += np.trapezoid(np.abs(Uref[j])**2,xref)
    field=np.sqrt(num/den)
    return float(qoi),float(field)


def split_marked(x, marked):
    marked=set(int(i) for i in marked)
    out=[x[0]]
    for e in range(len(x)-1):
        if e in marked: out.append(.5*(x[e]+x[e+1]))
        out.append(x[e+1])
    return np.array(out)


def build_adaptive(target_nodes:int, strategy:str, freqs:np.ndarray):
    x=np.linspace(0,L,13)
    while len(x)<target_nodes:
        m=assemble(x); U=solve_harmonic(m,freqs)
        if strategy=='snapshot':
            j=int(np.argmax(np.abs(U[:,-1])))
            eta=recovered_stress_indicator(m,U[j:j+1],freqs[j:j+1],np.ones(1))
        elif strategy=='pulse':
            tip=np.abs(U[:,-1])
            weights=.25+.75*tip/(tip.max()+1e-30)
            eta=recovered_stress_indicator(m,U,freqs,weights)
        else:
            raise ValueError(strategy)
        nadd=min(target_nodes-len(x),max(1,int(np.ceil(.30*(len(x)-1)))))
        marked=np.argsort(eta)[-nadd:]
        x=split_marked(x,marked)
    return x

FREQS=np.geomspace(8e3,280e3,30)

if __name__=='__main__':
    xref=np.linspace(0,L,2049)
    mref=assemble(xref)
    Uref=solve_harmonic(mref,FREQS)

    budgets=[13,21,33,49,73,97]
    rows=[]; meshes={}
    for n in budgets:
        for strategy in ('uniform','snapshot','pulse'):
            x=np.linspace(0,L,n) if strategy=='uniform' else build_adaptive(n,strategy,FREQS)
            m=assemble(x); U=solve_harmonic(m,FREQS)
            q,f=errors(m,U,xref,Uref,FREQS)
            rows.append((n,strategy,q,f,float(np.min(np.diff(x))),float(np.max(np.diff(x)))))
            meshes[(n,strategy)]=x

    print('nodes,strategy,qoi_rel_error,field_rel_error,h_min_mm,h_max_mm')
    for r in rows:
        print(f'{r[0]},{r[1]},{r[2]:.8f},{r[3]:.8f},{r[4]*1e3:.6f},{r[5]*1e3:.6f}')

    x=np.linspace(0,L,25); m=assemble(x)
    K=m.K[1:,1:]; F=m.F[1:]
    u=np.zeros(len(x)); u[1:]=spsolve(K,F)
    eta_snap=static_indicator(m,u)
    steps=64; amp=.35; eta_acc=np.zeros_like(eta_snap)
    for s in range(1,steps+1):
        scale=(s/steps)*(1+amp*np.sin(2*np.pi*s/steps))
        eta_acc += static_indicator(m,scale*u)**2
    eta_acc=np.sqrt(eta_acc)
    corr=np.corrcoef(eta_snap,eta_acc)[0,1]
    rank_equal=np.array_equal(np.argsort(eta_snap),np.argsort(eta_acc))
    print(f'quasistatic_indicator_correlation={corr:.12f}')
    print(f'quasistatic_rank_identical={rank_equal}')

    import json
    final_x=meshes[(49,'pulse')]
    final_m=assemble(final_x); final_u=solve_harmonic(final_m,FREQS)
    weights=.25+.75*np.abs(final_u[:,-1])/(np.abs(final_u[:,-1]).max()+1e-30)
    eta=recovered_stress_indicator(final_m,final_u,FREQS,weights)
    eta=(eta/(eta.max()+1e-30)).tolist()
    payload={
      'budgets':budgets,
      'rows':[{'nodes':n,'strategy':s,'qoiError':q,'fieldError':f,'hMinMm':h0*1e3,'hMaxMm':h1*1e3} for n,s,q,f,h0,h1 in rows],
      'quasiStaticCorrelation':corr,
      'quasiStaticRankIdentical':bool(rank_equal),
      'pulseMap49':{'nodes':final_x.tolist(),'indicator':eta},
      'frequenciesHz':FREQS.tolist()
    }
    with open('/mnt/data/piezo_sparse_benchmark_results.json','w') as fp: json.dump(payload,fp,indent=2)
