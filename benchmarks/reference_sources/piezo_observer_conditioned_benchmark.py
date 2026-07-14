from __future__ import annotations
import json, math, csv
from pathlib import Path
import numpy as np

DATA=json.loads(Path('/mnt/data/piezo_3d_nominal_data.json').read_text())
coords=np.asarray(DATA['coords_mm'],float)*1e-3
freqs=np.asarray(DATA['freqs'],float)
Uf=np.asarray(DATA['disp_full'],float); Ur=np.asarray(DATA['disp_reduced'],float)
Uf=Uf[...,0]+1j*Uf[...,1];Ur=Ur[...,0]+1j*Ur[...,1]
# shape freq,node,3
EPS0=8.8541878128e-12; NLAY=100; D33=400e-12;D31=-175e-12;E33=72e9*D33;E31=62e9*D31
VBASE=120.0

def rel(a,b):return float(np.linalg.norm(a-b)/(np.linalg.norm(b)+1e-30))
def phase_rms(a,b):
    wt=np.abs(b)/(np.max(np.abs(b))+1e-30);ph=np.angle(a*np.conjugate(b))
    return float(np.sqrt(np.sum(wt*ph*ph)/(np.sum(wt)+1e-30))*180/np.pi)
def worst(a,b,frac=.03):
    floor=frac*np.max(np.abs(b))+1e-30
    return float(np.max(np.abs(a-b)/(np.abs(b)+floor)))

pzt=[e for e in DATA['elements'] if e['material']=='PZT']
def avg_pzt_strains(U):
    out=[]
    for u in U:
        acc=np.zeros(3,complex);vt=0.0
        for e in pzt:
            nd=np.asarray(e['nodes'],int);c=coords[nd];du=u[nd]
            dx=c[:,0].max()-c[:,0].min();dy=c[:,1].max()-c[:,1].min();dz=c[:,2].max()-c[:,2].min();vol=dx*dy*dz
            ex=(du[[1,2,5,6],0].mean()-du[[0,3,4,7],0].mean())/(dx+1e-30)
            ey=(du[[2,3,6,7],1].mean()-du[[0,1,4,5],1].mean())/(dy+1e-30)
            ez=(du[[4,5,6,7],2].mean()-du[[0,1,2,3],2].mean())/(dz+1e-30)
            acc+=vol*np.array([ex,ey,ez]);vt+=vol
        out.append(acc/(vt+1e-30))
    return np.asarray(out)
SF=avg_pzt_strains(Uf);SR=avg_pzt_strains(Ur)
geom=DATA['geometry']; area=(2*geom['half_stack_mm']*1e-3)**2;h=geom['h_stack_mm']*1e-3

def electrical(epsr):
    C0=NLAY**2*EPS0*epsr*area/h
    qf=NLAY*area*(E31*(SF[:,0]+SF[:,1])+E33*SF[:,2])
    qr=NLAY*area*(E31*(SR[:,0]+SR[:,1])+E33*SR[:,2])
    w=2*np.pi*freqs
    Yf=1j*w*(C0+qf/VBASE);Yr=1j*w*(C0+qr/VBASE)
    return {'C0':C0,'qf':qf,'qr':qr,'Yf':Yf,'Yr':Yr,
            'admittance_rel':rel(Yr,Yf),'admittance_phase_deg':phase_rms(Yr,Yf),'admittance_worst':worst(Yr,Yf,.04),'motional_charge_rel':rel(qr,qf)}

top=np.where(np.isclose(coords[:,2],coords[:,2].max()))[0];xy=coords[top,:2]
def optical_fields(U,w0):
    S=U[:,top,2];w=np.exp(-2*((xy[:,0]/w0)**2+(xy[:,1]/w0)**2));w/=w.sum()+1e-30
    A=np.c_[np.ones(len(xy)),xy[:,0],xy[:,1]];W=np.sqrt(w)[:,None];Aw=W*A
    piston=[];tilt=[];wrms=[]
    for s in S:
        cr=np.linalg.lstsq(Aw,W[:,0]*s.real,rcond=None)[0];ci=np.linalg.lstsq(Aw,W[:,0]*s.imag,rcond=None)[0];c=cr+1j*ci
        res=s-A@c
        piston.append(c[0]);tilt.append(np.sqrt(np.abs(c[1])**2+np.abs(c[2])**2));wrms.append(np.sqrt(np.sum(w*np.abs(res)**2)))
    return np.asarray(piston),np.asarray(tilt),np.asarray(wrms)

def cavity_h1(piston,lamb,R1,R2,detuning,voltage):
    scale=voltage/VBASE;A=np.abs(piston)*scale;psi=np.angle(piston);beta=4*np.pi*A/lamb
    t=np.linspace(0,2*np.pi,96,endpoint=False);out=[];rr=np.sqrt(R1*R2);num=(1-R1)*(1-R2)
    for be,ps in zip(beta,psi):
        ph=detuning+be*np.cos(t+ps);power=num/(1+R1*R2-2*rr*np.cos(ph)+1e-30)
        out.append(2*np.mean(power*np.exp(-1j*t)))
    return np.asarray(out)

def optical(w0,lamb,R1,R2,detuning,voltage):
    pf,tf,wf=optical_fields(Uf,w0);pr,tr,wr=optical_fields(Ur,w0)
    scale=voltage/VBASE
    ovf=np.exp(-.5*((2*np.pi*w0/lamb)*tf*scale)**2-.5*(4*np.pi*wf*scale/lamb)**2)
    ovr=np.exp(-.5*((2*np.pi*w0/lamb)*tr*scale)**2-.5*(4*np.pi*wr*scale/lamb)**2)
    cf=cavity_h1(pf,lamb,R1,R2,detuning,voltage);cr=cavity_h1(pr,lamb,R1,R2,detuning,voltage)
    # wavefront denominator includes diffraction-relevant floor of lambda/10000
    wf_rel=float(np.linalg.norm(wr-wf)/(np.linalg.norm(wf)+np.sqrt(len(wf))*lamb/10000))
    return {'piston_rel':rel(pr,pf),'piston_phase_deg':phase_rms(pr,pf),'piston_worst':worst(pr,pf,.04),
            'wavefront_rel':wf_rel,'overlap_abs':float(np.max(np.abs(ovr-ovf))),
            'cavity_rel':rel(cr,cf),'cavity_phase_deg':phase_rms(cr,cf),'cavity_worst':worst(cr,cf,.05),
            'pf':pf,'pr':pr,'wf':wf,'wr':wr,'ovf':ovf,'ovr':ovr,'cf':cf,'cr':cr}

def certificate(e,o):
    E=max(e['admittance_rel']/.006,e['admittance_phase_deg']/.35,e['admittance_worst']/.025,e['motional_charge_rel']/.025)
    I=max(DATA['metrics']['port_rel']/.025,DATA['metrics']['port_shape']/.05,DATA['metrics']['bond_energy']/.08)
    C=max(o['piston_rel']/.005,o['piston_phase_deg']/.25,o['piston_worst']/.02,o['wavefront_rel']/.025,o['overlap_abs']/.015,o['cavity_rel']/.02,o['cavity_phase_deg']/.5,o['cavity_worst']/.05)
    return max(E,I,C),{'electromechanical_energy':E,'interface_information':I,'optical_coherence':C}

def main():
    epsrs=[800,1200,1600,2200]
    waists=[1.5e-3,2.5e-3,4e-3]
    lambdas=[532e-9,780e-9,1064e-9,1550e-9]
    voltages=[.25,1,5,20,80,120]
    cavities=[(.95,.999,.015),(.985,.9995,.006),(.995,.9998,.0025)]
    rows=[]
    e_cache={er:electrical(er) for er in epsrs}
    for er in epsrs:
      e=e_cache[er]
      for w0 in waists:
       for la in lambdas:
        for V in voltages:
         for R1,R2,det in cavities:
          o=optical(w0,la,R1,R2,det,V);delta,ch=certificate(e,o)
          rows.append({'epsr':er,'beam_waist_mm':w0*1e3,'wavelength_nm':la*1e9,'voltage':V,'R1':R1,'R2':R2,'detuning_rad':det,'delta':delta,'pass':delta<=1,**ch,
                       'admittance_rel':e['admittance_rel'],'admittance_phase_deg':e['admittance_phase_deg'],'admittance_worst':e['admittance_worst'],'motional_charge_rel':e['motional_charge_rel'],
                       **{k:o[k] for k in ['piston_rel','piston_phase_deg','piston_worst','wavefront_rel','overlap_abs','cavity_rel','cavity_phase_deg','cavity_worst']}})
    passes=sum(r['pass'] for r in rows);accepted=[r for r in rows if r['pass']]
    summary={'operating_points':len(rows),'reduced_passes':passes,'reduced_acceptance':passes/len(rows),'fallback_points':len(rows)-passes,
             'base_reduction':DATA['reduction'],'effective_median_exposed_fraction':(passes/len(rows))*(1-DATA['reduction'])+(1-passes/len(rows))*1,
             'max_delta':max(r['delta'] for r in rows),'median_delta':float(np.median([r['delta'] for r in rows])),
             'dominant_failures':{k:sum((not r['pass']) and r[k]>=max(r['electromechanical_energy'],r['interface_information'],r['optical_coherence'])-1e-12 for r in rows) for k in ['electromechanical_energy','interface_information','optical_coherence']},
             'emergent_statement':'Reduced-model validity is observer-conditioned. The physical state may be mechanically coherent yet require full-order fallback under readouts that amplify wavefront or cavity-sideband error.'}
    Path('/mnt/data/piezo_observer_conditioned_results.json').write_text(json.dumps({'summary':summary,'rows':rows},indent=2))
    with open('/mnt/data/piezo_observer_conditioned_results.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    # representative 1064 nm, 2.5 mm, 1V medium finesse for HTML
    e=electrical(1600);o=optical(2.5e-3,1064e-9,.985,.9995,.006,1);delta,ch=certificate(e,o)
    def cp(a):return [[float(z.real),float(z.imag)] for z in a]
    rep={'delta':delta,'channels':ch,'metrics':{k:v for k,v in {**e,**o}.items() if np.isscalar(v)},'admittance_full':cp(e['Yf']),'admittance_reduced':cp(e['Yr']),
         'cavity_full':cp(o['cf']),'cavity_reduced':cp(o['cr']),'piston_full':cp(o['pf']),'piston_reduced':cp(o['pr']),
         'wavefront_full':[float(x) for x in o['wf']],'wavefront_reduced':[float(x) for x in o['wr']],
         'overlap_full':[float(x) for x in o['ovf']],'overlap_reduced':[float(x) for x in o['ovr']],
         'observer':{'epsr':1600,'beam_waist_mm':2.5,'wavelength_nm':1064,'voltage':1,'R1':.985,'R2':.9995,'detuning_rad':.006},'sweep_summary':summary}
    Path('/mnt/data/piezo_observer_nominal.json').write_text(json.dumps(rep,separators=(',',':')))
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
