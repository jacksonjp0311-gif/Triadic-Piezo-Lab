from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
# import base benchmark module without running
import importlib.util,sys
spec=importlib.util.spec_from_file_location('ob','/mnt/data/piezo_observer_conditioned_benchmark.py')
ob=importlib.util.module_from_spec(spec);sys.modules['ob']=ob;spec.loader.exec_module(ob)

DATA=ob.DATA;freqs=ob.freqs;Uf=ob.Uf;Ur=ob.Ur;coords=ob.coords;top=ob.top;xy=ob.xy
# Precompute full/reduced top surface and electromechanical charge strain witness
Sf=Uf[:,top,2];Sr=Ur[:,top,2]
SF=ob.SF;SR=ob.SR
qf=ob.NLAY*ob.area*(ob.E31*(SF[:,0]+SF[:,1])+ob.E33*SF[:,2])
qr=ob.NLAY*ob.area*(ob.E31*(SR[:,0]+SR[:,1])+ob.E33*SR[:,2])
logf=np.log(freqs)

def interp_complex(vals, sent):
    sent=np.array(sorted(set(sent)),int);x=logf[sent]
    V=np.asarray(vals)
    shp=V.shape
    flat=V.reshape(len(V),-1)
    out=np.empty_like(flat,dtype=complex)
    for k in range(flat.shape[1]):
        out[:,k]=np.interp(logf,x,flat[sent,k].real)+1j*np.interp(logf,x,flat[sent,k].imag)
    return out.reshape(shp)

def optical_from_surface(S,w0):
    w=np.exp(-2*((xy[:,0]/w0)**2+(xy[:,1]/w0)**2));w/=w.sum()+1e-30
    A=np.c_[np.ones(len(xy)),xy[:,0],xy[:,1]];W=np.sqrt(w)[:,None];Aw=W*A
    piston=[];tilt=[];wrms=[]
    for s in S:
        cr=np.linalg.lstsq(Aw,W[:,0]*s.real,rcond=None)[0];ci=np.linalg.lstsq(Aw,W[:,0]*s.imag,rcond=None)[0];c=cr+1j*ci
        res=s-A@c;piston.append(c[0]);tilt.append(np.sqrt(np.abs(c[1])**2+np.abs(c[2])**2));wrms.append(np.sqrt(np.sum(w*np.abs(res)**2)))
    return np.asarray(piston),np.asarray(tilt),np.asarray(wrms)

def eval_all(sent):
    # Observer manifold stores only residual top-surface snapshots and terminal charge residual at sentinels.
    Sc=Sr+interp_complex(Sf-Sr,sent)
    qc=qr+interp_complex(qf-qr,sent)
    epsrs=[800,1200,1600,2200];waists=[1.5e-3,2.5e-3,4e-3];lambdas=[532e-9,780e-9,1064e-9,1550e-9];voltages=[.25,1,5,20,80,120];cavities=[(.95,.999,.015),(.985,.9995,.006),(.995,.9998,.0025)]
    rows=[];perfreq=np.zeros(len(freqs))
    for er in epsrs:
        C0=ob.NLAY**2*ob.EPS0*er*ob.area/ob.h;ww=2*np.pi*freqs
        Yf=1j*ww*(C0+qf/ob.VBASE);Yc=1j*ww*(C0+qc/ob.VBASE)
        e={'admittance_rel':ob.rel(Yc,Yf),'admittance_phase_deg':ob.phase_rms(Yc,Yf),'admittance_worst':ob.worst(Yc,Yf,.04),'motional_charge_rel':ob.rel(qc,qf)}
        ey=np.abs(Yc-Yf)/(np.abs(Yf)+.04*np.max(np.abs(Yf))+1e-30)/.025;perfreq=np.maximum(perfreq,ey)
        for w0 in waists:
            pf,tf,wf=optical_from_surface(Sf,w0);pc,tc,wc=optical_from_surface(Sc,w0)
            for la in lambdas:
                wf_rel=float(np.linalg.norm(wc-wf)/(np.linalg.norm(wf)+np.sqrt(len(wf))*la/10000))
                for V in voltages:
                    scale=V/ob.VBASE
                    ovf=np.exp(-.5*((2*np.pi*w0/la)*tf*scale)**2-.5*(4*np.pi*wf*scale/la)**2)
                    ovc=np.exp(-.5*((2*np.pi*w0/la)*tc*scale)**2-.5*(4*np.pi*wc*scale/la)**2)
                    for R1,R2,det in cavities:
                        cf=ob.cavity_h1(pf,la,R1,R2,det,V);cc=ob.cavity_h1(pc,la,R1,R2,det,V)
                        o={'piston_rel':ob.rel(pc,pf),'piston_phase_deg':ob.phase_rms(pc,pf),'piston_worst':ob.worst(pc,pf,.04),'wavefront_rel':wf_rel,'overlap_abs':float(np.max(np.abs(ovc-ovf))),'cavity_rel':ob.rel(cc,cf),'cavity_phase_deg':ob.phase_rms(cc,cf),'cavity_worst':ob.worst(cc,cf,.05)}
                        E=max(e['admittance_rel']/.006,e['admittance_phase_deg']/.35,e['admittance_worst']/.025,e['motional_charge_rel']/.025)
                        I=max(DATA['metrics']['port_rel']/.025,DATA['metrics']['port_shape']/.05,DATA['metrics']['bond_energy']/.08)
                        C=max(o['piston_rel']/.005,o['piston_phase_deg']/.25,o['piston_worst']/.02,o['wavefront_rel']/.025,o['overlap_abs']/.015,o['cavity_rel']/.02,o['cavity_phase_deg']/.5,o['cavity_worst']/.05)
                        delta=max(E,I,C);rows.append(delta)
                        cv=np.abs(cc-cf)/(np.abs(cf)+.05*np.max(np.abs(cf))+1e-30)/.05
                        pp=np.abs(pc-pf)/(np.abs(pf)+.04*np.max(np.abs(pf))+1e-30)/.02
                        perfreq=np.maximum(perfreq,np.maximum(cv,pp))
    return {'acceptance':float(np.mean(np.array(rows)<=1)),'max_delta':float(np.max(rows)),'median_delta':float(np.median(rows)),'perfreq':perfreq,'Sc':Sc,'qc':qc}

def main():
    sent=sorted(set(DATA['sentinels']));history=[]
    for it in range(24):
        r=eval_all(sent);history.append({'iteration':it+1,'sentinels':sent.copy(),'count':len(sent),'acceptance':r['acceptance'],'max_delta':r['max_delta'],'median_delta':r['median_delta']})
        print(it+1,len(sent),round(r['acceptance'],4),round(r['median_delta'],3),round(r['max_delta'],2),flush=True)
        if r['acceptance']>=.995 or len(sent)>=24:break
        candidates=[i for i in range(len(freqs)) if i not in sent]
        worst=max(candidates,key=lambda i:r['perfreq'][i]);sent=sorted(sent+[worst])
    final=eval_all(sent)
    # compact representative corrected witness at 1064nm / 2.5mm / 1V
    pf,tf,wf=optical_from_surface(Sf,2.5e-3);pc,tc,wc=optical_from_surface(final['Sc'],2.5e-3)
    cf=ob.cavity_h1(pf,1064e-9,.985,.9995,.006,1);cc=ob.cavity_h1(pc,1064e-9,.985,.9995,.006,1)
    er=1600;C0=ob.NLAY**2*ob.EPS0*er*ob.area/ob.h;ww=2*np.pi*freqs;Yf=1j*ww*(C0+qf/ob.VBASE);Yc=1j*ww*(C0+final['qc']/ob.VBASE)
    def cp(a):return [[float(z.real),float(z.imag)] for z in a]
    payload={'method':'observer-manifold adversarial closure','base_global_dofs':DATA['global_dofs'],'full_dofs':DATA['full_dofs'],'base_reduction':DATA['reduction'],'observer_surface_coordinates':len(top),'final_sentinels':sent,'history':history,
             'final':{'acceptance':final['acceptance'],'max_delta':final['max_delta'],'median_delta':final['median_delta'],'observer_snapshot_scalars':len(sent)*(2*len(top)+2),
                      'emergent_statement':'When a nonlinear observer amplifies small state errors, stability can be restored by a sparse calibrated witness manifold without reopening the hidden volume. The observer manifold must itself be adversarially closed over frequency.'},
             'representative':{'admittance_full':cp(Yf),'admittance_corrected':cp(Yc),'cavity_full':cp(cf),'cavity_corrected':cp(cc),'piston_full':cp(pf),'piston_corrected':cp(pc),'wavefront_full':[float(x) for x in wf],'wavefront_corrected':[float(x) for x in wc]}}
    Path('/mnt/data/piezo_observer_manifold_results.json').write_text(json.dumps(payload,indent=2))
    print(json.dumps(payload['final'],indent=2))
if __name__=='__main__':main()
