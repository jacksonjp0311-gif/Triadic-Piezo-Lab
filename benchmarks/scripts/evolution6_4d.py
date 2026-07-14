#!/usr/bin/env python3
from __future__ import annotations
import json, math, csv
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from numba import njit
OUT=Path(__import__('os').environ.get('TRIAD_BENCH_OUT', Path(__file__).resolve().parents[1]/'results'/'evolution6_4d'/'reproduced'));OUT.mkdir(parents=True,exist_ok=True)

@njit(cache=True)
def floquet_exp(zeta,h,ratio,steps=220):
    T=math.pi/ratio;dt=T/steps
    y00=1.;y01=0.;y10=0.;y11=1.;t=0.
    for _ in range(steps):
        # AY helper expanded for each stage
        k=1+h*math.cos(2*ratio*t)
        a00=y10;a01=y11;a10=-k*y00-2*zeta*y10;a11=-k*y01-2*zeta*y11
        t2=t+dt*.5
        z00=y00+dt*a00*.5;z01=y01+dt*a01*.5;z10=y10+dt*a10*.5;z11=y11+dt*a11*.5
        k=1+h*math.cos(2*ratio*t2)
        b00=z10;b01=z11;b10=-k*z00-2*zeta*z10;b11=-k*z01-2*zeta*z11
        z00=y00+dt*b00*.5;z01=y01+dt*b01*.5;z10=y10+dt*b10*.5;z11=y11+dt*b11*.5
        c00=z10;c01=z11;c10=-k*z00-2*zeta*z10;c11=-k*z01-2*zeta*z11
        t4=t+dt
        z00=y00+dt*c00;z01=y01+dt*c01;z10=y10+dt*c10;z11=y11+dt*c11
        k=1+h*math.cos(2*ratio*t4)
        d00=z10;d01=z11;d10=-k*z00-2*zeta*z10;d11=-k*z01-2*zeta*z11
        y00+=dt*(a00+2*b00+2*c00+d00)/6;y01+=dt*(a01+2*b01+2*c01+d01)/6
        y10+=dt*(a10+2*b10+2*c10+d10)/6;y11+=dt*(a11+2*b11+2*c11+d11)/6;t=t4
    tr=y00+y11;det=y00*y11-y01*y10;disc=tr*tr-4*det
    if disc>=0:
        root=math.sqrt(disc);l1=.5*(tr+root);l2=.5*(tr-root);rho=max(abs(l1),abs(l2))
    else:rho=math.sqrt(abs(det))
    return math.log(max(rho,1e-300))/T

def floquet_threshold(zeta,ratio=1.):
    f=lambda h:float(floquet_exp(zeta,h,ratio))
    lo=0.;hi=max(.08,10*zeta)
    while f(hi)<=0 and hi<2:hi*=1.7
    return float(brentq(f,lo,hi,xtol=1e-10))

@njit(cache=True)
def duffing_sweep_numba(freqs,zeta,beta,force,reverse,periods=150,steps=105):
    nfreq=len(freqs);amps=np.zeros(nfreq);phases=np.zeros(nfreq);x=0.;v=0.;phase=0.
    for kk in range(nfreq):
        idx=nfreq-1-kk if reverse else kk;om=freqs[idx];T=2*math.pi/om;dt=T/steps;keep_periods=30
        sc=0.;ss=0.;cnt=0
        for n in range(periods*steps):
            # phase-continuous forcing across the frequency sweep
            p1=phase
            k1x=v;k1v=force*math.cos(p1)-2*zeta*v-x-beta*x*x*x
            x2=x+dt*k1x*.5;v2=v+dt*k1v*.5;p2=phase+om*dt*.5
            k2x=v2;k2v=force*math.cos(p2)-2*zeta*v2-x2-beta*x2*x2*x2
            x3=x+dt*k2x*.5;v3=v+dt*k2v*.5
            k3x=v3;k3v=force*math.cos(p2)-2*zeta*v3-x3-beta*x3*x3*x3
            x4=x+dt*k3x;v4=v+dt*k3v;p4=phase+om*dt
            k4x=v4;k4v=force*math.cos(p4)-2*zeta*v4-x4-beta*x4*x4*x4
            x+=dt*(k1x+2*k2x+2*k3x+k4x)/6;v+=dt*(k1v+2*k2v+2*k3v+k4v)/6;phase=p4
            if n>=(periods-keep_periods)*steps:
                sc+=x*math.cos(phase);ss+=x*math.sin(phase);cnt+=1
        c=2*sc/cnt;s=2*ss/cnt;amps[idx]=math.sqrt(c*c+s*s);phases[idx]=math.atan2(-s,c)
    return amps,phases

def duffing_roots(omega,zeta,beta,force):
    d=1-omega*omega;q=.75*beta;c=2*zeta*omega
    roots=np.roots([q*q,2*d*q,d*d+c*c,-force*force]);ans=[]
    for r in roots:
        if abs(r.imag)<1e-8 and r.real>0:
            y=float(r.real);der=3*q*q*y*y+4*d*q*y+(d*d+c*c);ans.append((math.sqrt(y),der>0))
    return sorted(ans)
def hb_sweep(freqs,zeta,beta,force,reverse=False):
    order=range(len(freqs)-1,-1,-1) if reverse else range(len(freqs));out=np.zeros(len(freqs));prev=0.
    for idx in order:
        allr=duffing_roots(float(freqs[idx]),zeta,beta,force);cand=[a for a,s in allr if s] or [a for a,_ in allr]
        chosen=(max(cand) if reverse else min(cand)) if prev==0 else min(cand,key=lambda a:abs(a-prev));out[idx]=chosen;prev=chosen
    return out

@njit(cache=True)
def cavity_core(ratio,A=.22,g=2.2,det=-.65,drive=1.,omega=1.,periods=65,steps=120):
    kappa=ratio*omega;T=2*math.pi/omega;dt=T/steps;ar=0.;ai=0.;t=0.;keep=22*steps
    Id=np.zeros(keep);Iad=np.zeros(keep);Ico=np.zeros(keep);tt=np.zeros(keep);j=0
    total=periods*steps
    for n in range(total):
        # complex derivative expanded: da=(iD-k/2)a+drive
        x=A*math.cos(omega*t);D=det+g*x
        k1r=-kappa*.5*ar-D*ai+drive;k1i=D*ar-kappa*.5*ai
        t2=t+dt*.5;x=A*math.cos(omega*t2);D=det+g*x
        ar2=ar+dt*k1r*.5;ai2=ai+dt*k1i*.5;k2r=-kappa*.5*ar2-D*ai2+drive;k2i=D*ar2-kappa*.5*ai2
        ar3=ar+dt*k2r*.5;ai3=ai+dt*k2i*.5;k3r=-kappa*.5*ar3-D*ai3+drive;k3i=D*ar3-kappa*.5*ai3
        t4=t+dt;x=A*math.cos(omega*t4);D=det+g*x;ar4=ar+dt*k3r;ai4=ai+dt*k3i
        k4r=-kappa*.5*ar4-D*ai4+drive;k4i=D*ar4-kappa*.5*ai4
        ar+=dt*(k1r+2*k2r+2*k3r+k4r)/6;ai+=dt*(k1i+2*k2i+2*k3i+k4i)/6;t=t4
        if n>=total-keep:
            x=A*math.cos(omega*t);xd=-A*omega*math.sin(omega*t);D=det+g*x;lr=kappa*.5;li=-D;den=lr*lr+li*li
            a0r=drive*lr/den;a0i=-drive*li/den
            # da_ss/dt = i*g*xd*drive/L^2; complex arithmetic
            # 1/L^2 = ((lr-li i)^2)/den^2
            qrr=(lr*lr-li*li)/(den*den);qii=-2*lr*li/(den*den)
            dr=-g*xd*drive*qii;di=g*xd*drive*qrr
            # correction = derivative/L
            cr=(dr*lr+di*li)/den;ci=(di*lr-dr*li)/den
            acr=a0r-cr;aci=a0i-ci
            Id[j]=ar*ar+ai*ai;Iad[j]=a0r*a0r+a0i*a0i;Ico[j]=acr*acr+aci*aci;tt[j]=t;j+=1
    return Id,Iad,Ico,tt

def harmonic(y,t,omega):return 2/len(y)*np.sum(y*np.exp(-1j*omega*t))
def cavity_metrics(r):
    Id,Ia,Ic,t=cavity_core(r);Hd=harmonic(Id,t,1.);Ha=harmonic(Ia,t,1.);Hc=harmonic(Ic,t,1.)
    def m(H):return float(abs(H-Hd)/(abs(Hd)+1e-14)),float(abs(np.angle(H/Hd,deg=True)))
    ea,pa=m(Ha);ec,pc=m(Hc)
    return dict(kappa_over_omega=r,adiabatic_harmonic_error=ea,adiabatic_phase_error_deg=pa,corrected_harmonic_error=ec,corrected_phase_error_deg=pc,
                adiabatic_waveform_error=float(np.linalg.norm(Ia-Id)/np.linalg.norm(Id)),corrected_waveform_error=float(np.linalg.norm(Ic-Id)/np.linalg.norm(Id)))

@njit(cache=True)
def coupled_core(mode,ratio,zeta=.022,beta=.15,force=.055,omega=1.03,Aopt=.035,g=1.8,det=-.7,drive=1.,periods=135,steps=95):
    kappa=ratio*omega;T=2*math.pi/omega;dt=T/steps;x=0.;v=0.;ar=0.;ai=0.;t=0.;keep=28*steps
    xs=np.zeros(keep);Is=np.zeros(keep);ts=np.zeros(keep);j=0;total=periods*steps
    for n in range(total):
        # derivative function manually, mode 0 dynamic 1 adiabatic 2 corrected
        vals=np.zeros(16)
        for stage in range(4):
            if stage==0:tt=t;xx=x;vv=v;aar=ar;aai=ai
            elif stage==1:tt=t+dt*.5;xx=x+dt*vals[0]*.5;vv=v+dt*vals[1]*.5;aar=ar+dt*vals[2]*.5;aai=ai+dt*vals[3]*.5
            elif stage==2:tt=t+dt*.5;xx=x+dt*vals[4]*.5;vv=v+dt*vals[5]*.5;aar=ar+dt*vals[6]*.5;aai=ai+dt*vals[7]*.5
            else:tt=t+dt;xx=x+dt*vals[8];vv=v+dt*vals[9];aar=ar+dt*vals[10];aai=ai+dt*vals[11]
            D=det+g*xx
            if mode==0:or_=aar;oi=aai;dar=-kappa*.5*aar-D*aai+drive;dai=D*aar-kappa*.5*aai
            else:
                lr=kappa*.5;li=-D;den=lr*lr+li*li;or_=drive*lr/den;oi=-drive*li/den
                if mode==2:
                    qrr=(lr*lr-li*li)/(den*den);qii=-2*lr*li/(den*den);dr=-g*vv*drive*qii;di=g*vv*drive*qrr
                    cr=(dr*lr+di*li)/den;ci=(di*lr-dr*li)/den;or_-=cr;oi-=ci
                dar=0.;dai=0.
            dx=vv;dv=force*math.cos(omega*tt)-2*zeta*vv-xx-beta*xx*xx*xx+Aopt*(or_*or_+oi*oi)
            base=stage*4;vals[base]=dx;vals[base+1]=dv;vals[base+2]=dar;vals[base+3]=dai
        x+=dt*(vals[0]+2*vals[4]+2*vals[8]+vals[12])/6;v+=dt*(vals[1]+2*vals[5]+2*vals[9]+vals[13])/6
        if mode==0:
            ar+=dt*(vals[2]+2*vals[6]+2*vals[10]+vals[14])/6;ai+=dt*(vals[3]+2*vals[7]+2*vals[11]+vals[15])/6
        t+=dt
        if n>=total-keep:
            D=det+g*x
            if mode==0:or_=ar;oi=ai
            else:
                lr=kappa*.5;li=-D;den=lr*lr+li*li;or_=drive*lr/den;oi=-drive*li/den
                if mode==2:
                    qrr=(lr*lr-li*li)/(den*den);qii=-2*lr*li/(den*den);dr=-g*v*drive*qii;di=g*v*drive*qrr
                    cr=(dr*lr+di*li)/den;ci=(di*lr-dr*li)/den;or_-=cr;oi-=ci
            xs[j]=x;Is[j]=or_*or_+oi*oi;ts[j]=t;j+=1
    return xs,Is,ts

def relphase(H,ref):return float(abs(H-ref)/(abs(ref)+1e-14)),float(abs(np.angle(H/ref,deg=True)))

def main():
    # warm compile
    floquet_exp(.02,.1,1.);duffing_sweep_numba(np.array([1.]),.02,.1,.03,False);cavity_core(2.);coupled_core(0,2.)
    zetas=[.005,.01,.02,.04,.07];thresholds=[]
    for z in zetas:
        hc=floquet_threshold(z);thresholds.append(dict(zeta=z,numerical_hcrit=hc,asymptotic_4zeta=4*z,relative_error=abs(hc-4*z)/hc))
    ratios=np.linspace(.82,1.18,35);hs=np.linspace(0,.34,35);fmap=np.empty((len(hs),len(ratios)))
    for i,h in enumerate(hs):
        for j,r in enumerate(ratios):fmap[i,j]=floquet_exp(.02,float(h),float(r),180)
    # low-amplitude validation where first-harmonic balance should be single-valued
    freqs_val=np.linspace(.82,1.18,33);zeta_val=.03;beta_val=.06;force_val=.018
    dv,_=duffing_sweep_numba(freqs_val,zeta_val,beta_val,force_val,False);hv=hb_sweep(freqs_val,zeta_val,beta_val,force_val)
    eval_=float(np.linalg.norm(hv-dv)/np.linalg.norm(dv))
    # stronger nonlinear sweep to reveal branch history and jump behavior
    freqs=np.linspace(.92,1.13,57);zeta=.025;beta=.22;force=.065
    du,_=duffing_sweep_numba(freqs,zeta,beta,force,False);dd,_=duffing_sweep_numba(freqs,zeta,beta,force,True);hu=hb_sweep(freqs,zeta,beta,force);hd=hb_sweep(freqs,zeta,beta,force,True)
    eu=float(np.linalg.norm(hu-du)/np.linalg.norm(du));ed=float(np.linalg.norm(hd-dd)/np.linalg.norm(dd));area_d=float(np.trapezoid(abs(du-dd),freqs));area_h=float(np.trapezoid(abs(hu-hd),freqs))
    cav=[cavity_metrics(r) for r in [.35,.5,.75,1,1.5,2,3,5,8,12,20]]
    coupled=[];sample=None
    for r in [.5,1,2,4,8,16]:
        xd,Id,t=coupled_core(0,r);xa,Ia,_=coupled_core(1,r);xc,Ic,_=coupled_core(2,r)
        Hxd=harmonic(xd,t,1.03);Hid=harmonic(Id,t,1.03);Hxa=harmonic(xa,t,1.03);Hia=harmonic(Ia,t,1.03);Hxc=harmonic(xc,t,1.03);Hic=harmonic(Ic,t,1.03)
        exa,pxa=relphase(Hxa,Hxd);eia,pia=relphase(Hia,Hid);exc,pxc=relphase(Hxc,Hxd);eic,pic=relphase(Hic,Hid)
        coupled.append(dict(kappa_over_omega=r,adiabatic_mech_error=exa,adiabatic_mech_phase_deg=pxa,adiabatic_optical_error=eia,adiabatic_optical_phase_deg=pia,corrected_mech_error=exc,corrected_mech_phase_deg=pxc,corrected_optical_error=eic,corrected_optical_phase_deg=pic))
        if r==2:sample={'t':t[-190:].tolist(),'x':xd[-190:].tolist(),'intensity':Id[-190:].tolist()}
    results={'model_scope':'dimensionless reduced-order validation; not device calibration','floquet_thresholds':thresholds,
      'floquet_map':{'ratios':ratios.tolist(),'h':hs.tolist(),'exponent':fmap.tolist(),'zeta':.02},
      'duffing':{'validation_frequency_ratio':freqs_val.tolist(),'validation_direct':dv.tolist(),'validation_hb':hv.tolist(),'validation_relative_error':eval_,'frequency_ratio':freqs.tolist(),'direct_up':du.tolist(),'direct_down':dd.tolist(),'hb_up':hu.tolist(),'hb_down':hd.tolist(),'relative_error_up':eu,'relative_error_down':ed,'hysteresis_area_direct':area_d,'hysteresis_area_hb':area_h,'zeta':zeta,'beta':beta,'force':force},
      'cavity_memory':cav,'coupled_4d':coupled,'trajectory_sample':sample,
      'emergent':{'principle':'observer coherence requires a memory state when cavity bandwidth is not large compared with mechanical frequency','stability_gate':'Floquet growth, nonlinear branch ambiguity, and optical-memory mismatch must be certified independently'}}
    (OUT/'piezo_4d_floquet_results.json').write_text(json.dumps(results,indent=2))
    rows=[]
    for x in thresholds:rows.append({'benchmark':'floquet_threshold',**x})
    for x in cav:rows.append({'benchmark':'cavity_memory',**x})
    for x in coupled:rows.append({'benchmark':'coupled_4d',**x})
    keys=sorted({k for r in rows for k in r})
    with (OUT/'piezo_4d_floquet_results.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
    summary={'floquet_max_relative_threshold_error':max(x['relative_error'] for x in thresholds),'duffing_validation_relative_error':eval_,'duffing_relative_error_up':eu,'duffing_relative_error_down':ed,'duffing_hysteresis_area_ratio':area_h/area_d,'adiabatic_pass_count_5pct':sum(x['adiabatic_harmonic_error']<.05 for x in cav),'corrected_pass_count_5pct':sum(x['corrected_harmonic_error']<.05 for x in cav),'coupled_adiabatic_optical_pass_count_5pct':sum(x['adiabatic_optical_error']<.05 for x in coupled),'coupled_corrected_optical_pass_count_5pct':sum(x['corrected_optical_error']<.05 for x in coupled)}
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
