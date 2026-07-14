from __future__ import annotations

import csv
import json
import math
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import spsolve

OUT = Path('/mnt/data')
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


@dataclass(frozen=True)
class Scenario:
    name: str
    freqs: np.ndarray
    e_profile: Callable[[np.ndarray], np.ndarray]
    a_profile: Callable[[np.ndarray], np.ndarray]
    description: str


@dataclass
class Model:
    x: np.ndarray
    K: csr_matrix
    M: csr_matrix
    C: csr_matrix
    F: np.ndarray
    scenario: Scenario


@dataclass
class Result:
    cycle: int
    scenario: str
    nodes: int
    elements: int
    strategy: str
    qoi_rel_error: float
    field_rel_error: float
    peak_frequency_rel_error: float
    h_min_mm: float
    h_max_mm: float
    h_ratio: float
    locality_index: float
    coverage_demand: float
    build_runtime_ms: float


def gaussian(x, c, s):
    return np.exp(-((x-c)/s)**2)


def make_profiles(kind: str):
    if kind == 'standard':
        def ef(x):
            x=np.asarray(x); return E0*np.clip(1-.58*gaussian(x,.0142,.00062)-.28*gaussian(x,.0180,.00038)-.055*np.sin(2*np.pi*x/.0025)**2,.24,None)
        def af(x):
            x=np.asarray(x); return A0*(1-.42*gaussian(x,.0185,.00058))
    elif kind == 'shifted':
        def ef(x):
            x=np.asarray(x); return E0*np.clip(1-.52*gaussian(x,.0068,.00052)-.32*gaussian(x,.0104,.00044)-.04*np.sin(2*np.pi*x/.0031)**2,.28,None)
        def af(x):
            x=np.asarray(x); return A0*(1-.36*gaussian(x,.0111,.00066))
    elif kind == 'dual':
        def ef(x):
            x=np.asarray(x); return E0*np.clip(1-.40*gaussian(x,.0052,.00042)-.48*gaussian(x,.0147,.00055)-.22*gaussian(x,.0182,.00034),.25,None)
        def af(x):
            x=np.asarray(x); return A0*(1-.22*gaussian(x,.0060,.00050)-.38*gaussian(x,.0185,.00060))
    elif kind == 'smooth':
        def ef(x):
            x=np.asarray(x); return E0*(.88+.12*np.cos(np.pi*x/L)**2)
        def af(x):
            x=np.asarray(x); return A0*(.90+.10*np.cos(2*np.pi*x/L)**2)
    else:
        raise ValueError(kind)
    return ef, af


def make_scenarios():
    scenarios=[]
    ef,af=make_profiles('standard')
    scenarios.append(Scenario('localized_standard_3_60kHz',np.geomspace(3e3,60e3,24),ef,af,'Compliant bond and mirror neck; sub-resonant/localized interrogation.'))
    ef,af=make_profiles('shifted')
    scenarios.append(Scenario('localized_shifted_4_80kHz',np.geomspace(4e3,80e3,24),ef,af,'Held-out shifted defect geometry tests spatial generalization.'))
    ef,af=make_profiles('dual')
    scenarios.append(Scenario('dual_interface_5_120kHz',np.geomspace(5e3,120e3,28),ef,af,'Two separated compliance zones test multi-hotspot allocation.'))
    ef,af=make_profiles('standard')
    scenarios.append(Scenario('broadband_standard_8_280kHz',np.geomspace(8e3,280e3,32),ef,af,'Broadband wave field where global resolution may be required.'))
    ef,af=make_profiles('smooth')
    scenarios.append(Scenario('smooth_global_12_300kHz',np.geomspace(12e3,300e3,32),ef,af,'No localized defect; negative control for inappropriate sparsening.'))
    return scenarios


def assemble(x: np.ndarray, scenario: Scenario) -> Model:
    n=len(x); K=lil_matrix((n,n),dtype=float); M=lil_matrix((n,n),dtype=float); F=np.zeros(n)
    eps=D33*VOLTAGE*N_LAYERS/L
    for e in range(n-1):
        a,b=x[e],x[e+1]; h=b-a; xm=.5*(a+b)
        E=float(scenario.e_profile(np.array([xm]))[0]); A=float(scenario.a_profile(np.array([xm]))[0]); rho=RHO0
        c=E*A/h; mc=rho*A*h/6; fe=E*A*eps
        i,j=e,e+1
        K[i,i]+=c; K[i,j]-=c; K[j,i]-=c; K[j,j]+=c
        M[i,i]+=2*mc; M[i,j]+=mc; M[j,i]+=mc; M[j,j]+=2*mc
        F[i]-=fe; F[j]+=fe
    M[-1,-1]+=MIRROR_MASS; K[-1,-1]+=MIRROR_K
    K=K.tocsr(); M=M.tocsr()
    w1=2*np.pi*28e3; w2=2*np.pi*210e3
    beta=2*ZETA/(w1+w2); alpha=beta*w1*w2
    C=(alpha*M+beta*K).tocsr()
    return Model(x,K,M,C,F,scenario)


def solve_harmonic(model: Model, freqs: np.ndarray) -> np.ndarray:
    K=model.K[1:,1:]; M=model.M[1:,1:]; C=model.C[1:,1:]; F=model.F[1:]
    U=np.zeros((len(freqs),len(model.x)),complex)
    for j,f in enumerate(freqs):
        w=2*np.pi*f; Z=(K-w*w*M+1j*w*C).tocsc(); U[j,1:]=spsolve(Z,F)
    return U


def element_indicator_matrix(model: Model, U: np.ndarray, freqs: np.ndarray) -> np.ndarray:
    x=model.x; h=np.diff(x); mid=.5*(x[:-1]+x[1:]); E=model.scenario.e_profile(mid); A=model.scenario.a_profile(mid); rho=RHO0
    out=np.zeros((len(freqs),len(h)))
    for j,f in enumerate(freqs):
        u=U[j]; sig=E*np.diff(u)/h
        sn=np.zeros(len(x),complex); sn[0]=sig[0]; sn[-1]=sig[-1]
        for i in range(1,len(x)-1):
            wl=h[i-1]; wr=h[i]; sn[i]=(wr*sig[i-1]+wl*sig[i])/(wl+wr)
        diffL=sn[:-1]-sig; diffR=sn[1:]-sig
        recovery=h*A/np.maximum(E,1.0)*.5*(np.abs(diffL)**2+np.abs(diffR)**2)
        um=.5*(u[:-1]+u[1:]); w=2*np.pi*f
        dyn=np.abs(rho*A*w*w*um); residual=h**3/np.maximum(E*A,1e-30)*dyn**2/12
        # Phase/curvature channel: catches propagating waves missed by purely local recovery.
        phase_curv=np.zeros_like(h)
        if len(u)>2:
            du=np.diff(u)/h
            jump=np.zeros_like(h,dtype=float)
            jump[:-1]+=np.abs(du[1:]-du[:-1])**2
            jump[1:]+=np.abs(du[1:]-du[:-1])**2
            phase_curv=(h**2)*jump
        out[j]=np.maximum(recovery+.12*residual+.018*phase_curv,0)
    return out


def normalize(v):
    v=np.asarray(v,float); mx=float(np.max(v)) if v.size else 0.0
    return v/(mx+1e-30)


def spatial_locality(energy: np.ndarray):
    q=np.maximum(np.asarray(energy,float),0); total=q.sum()
    if total<=1e-30: return 0.0
    p=q/total; ent=-np.sum(p*np.log(p+1e-30)); neff=float(np.exp(ent)); n=len(q)
    return float(np.clip(1-neff/max(n,1),0,1))


def pulse_data(model: Model, U: np.ndarray, freqs: np.ndarray):
    E2=element_indicator_matrix(model,U,freqs)
    tip=np.abs(U[:,-1]); weights=.30+.70*tip/(tip.max()+1e-30)
    agg=np.sum(weights[:,None]*E2,axis=0)
    eta=np.sqrt(np.maximum(agg,0))
    locality=spatial_locality(agg)
    # Spectral novelty: elements that become important only in part of the pulse.
    norm_rows=E2/(np.max(E2,axis=1,keepdims=True)+1e-30)
    novelty=np.std(norm_rows,axis=0)
    return eta,normalize(novelty),locality,E2,weights


def static_snapshot_data(model,U,freqs):
    j=int(np.argmax(np.abs(U[:,-1]))); E2=element_indicator_matrix(model,U[j:j+1],freqs[j:j+1])
    return np.sqrt(E2[0]),0.0


def wave_coverage(model: Model, freqs: np.ndarray, points_per_wavelength=8.0):
    x=model.x; h=np.diff(x); mid=.5*(x[:-1]+x[1:]); c=np.sqrt(model.scenario.e_profile(mid)/RHO0)
    h_allow=c/(float(np.max(freqs))*points_per_wavelength)
    pressure=np.maximum(h/(h_allow+1e-30)-1,0)
    demand=float(np.mean(pressure>0))
    return normalize(pressure),demand,h_allow


def neighbor_smooth(v):
    v=np.asarray(v,float); out=.58*v.copy()
    if len(v)>1:
        out[:-1]+=.21*v[1:]; out[1:]+=.21*v[:-1]
    return out


def split_marked(x, marked):
    marked=set(int(i) for i in marked); out=[x[0]]
    for e in range(len(x)-1):
        if e in marked: out.append(.5*(x[e]+x[e+1]))
        out.append(x[e+1])
    return np.array(out)


def choose_marked(score: np.ndarray, energy: np.ndarray, capacity: int, theta=.65):
    order=np.argsort(score)[::-1]; total=float(np.sum(np.maximum(energy,0))); picked=[]; acc=0.0
    for i in order:
        picked.append(int(i)); acc+=max(float(energy[i]),0.0)
        if len(picked)>=capacity: break
        if total>0 and acc>=theta*total and len(picked)>=1: break
    return np.array(picked,dtype=int)


def build_adaptive(target_nodes: int, strategy: str, scenario: Scenario, cycle: int):
    x=np.linspace(0,L,13); locality=0.0; demand=0.0
    while len(x)<target_nodes:
        model=assemble(x,scenario); U=solve_harmonic(model,scenario.freqs); h=np.diff(x)
        capacity=min(target_nodes-len(x),max(1,int(math.ceil(.22*len(h)))))
        wave,demand,_=wave_coverage(model,scenario.freqs,8.0)
        size=normalize(h)
        if strategy=='snapshot':
            eta,_=static_snapshot_data(model,U,scenario.freqs); score=normalize(eta); energy=eta**2
        elif strategy=='pulse_v0':
            eta,novelty,locality,E2,weights=pulse_data(model,U,scenario.freqs); score=normalize(eta); energy=eta**2
        elif strategy=='cycle1_regularized':
            eta,novelty,locality,E2,weights=pulse_data(model,U,scenario.freqs)
            local=neighbor_smooth(normalize(eta))
            # Cycle 1: prevent starvation with explicit wave and size pressure.
            score=.68*local+.20*wave+.12*size+.08*novelty
            energy=eta**2+.20*np.max(eta**2)*wave+.05*np.max(eta**2)*size
        elif strategy=='cycle2_governed':
            eta,novelty,locality,E2,weights=pulse_data(model,U,scenario.freqs)
            local=neighbor_smooth(normalize(eta))
            # Cycle 2: entropy-controlled dual channel. Local residual channel receives
            # authority only when the pulse map is concentrated. Coverage receives the
            # remaining budget and closes the compression gate for global fields.
            authority=float(np.clip((locality-.08)/.48,0,1))
            coverage=np.maximum(size,wave)
            score=authority*(.82*local+.18*novelty)+(1-authority)*coverage
            # Always retain a small counter-channel so neither local nor global evidence can vanish.
            score+=.10*local+.08*coverage
            energy=authority*(eta**2)+(1-authority)*(np.max(eta**2)+1e-30)*coverage**2
        else:
            raise ValueError(strategy)
        marked=choose_marked(score,energy,capacity,theta=.68 if cycle==1 else .62)
        x=split_marked(x,marked)
    return x,locality,demand


def interp_complex(xc,uc,xr):
    return np.interp(xr,xc,uc.real)+1j*np.interp(xr,xc,uc.imag)


def errors(model,U,xref,Uref,freqs):
    tip=U[:,-1]; tipref=Uref[:,-1]
    qoi=float(np.linalg.norm(tip-tipref)/(np.linalg.norm(tipref)+1e-30))
    num=0.0; den=0.0
    for j in range(len(freqs)):
        ui=interp_complex(model.x,U[j],xref)
        num+=np.trapezoid(np.abs(ui-Uref[j])**2,xref); den+=np.trapezoid(np.abs(Uref[j])**2,xref)
    field=float(np.sqrt(num/(den+1e-30)))
    peak=float(abs(freqs[int(np.argmax(np.abs(tip)))]-freqs[int(np.argmax(np.abs(tipref)))])/(freqs[int(np.argmax(np.abs(tipref)))]+1e-30))
    return qoi,field,peak


def run_cycle(cycle: int, strategies: Iterable[str], scenarios: list[Scenario], budgets=(21,33,49,73)):
    rows=[]; maps={}; start=time.perf_counter()
    for scenario in scenarios:
        xref=np.linspace(0,L,1025); mref=assemble(xref,scenario); Uref=solve_harmonic(mref,scenario.freqs)
        for n in budgets:
            for strategy in strategies:
                t0=time.perf_counter()
                if strategy=='uniform':
                    x=np.linspace(0,L,n); loc=0.0; demand=0.0
                else:
                    x,loc,demand=build_adaptive(n,strategy,scenario,cycle)
                model=assemble(x,scenario); U=solve_harmonic(model,scenario.freqs); q,f,p=errors(model,U,xref,Uref,scenario.freqs)
                h=np.diff(x)
                row=Result(cycle,scenario.name,n,n-1,strategy,q,f,p,float(h.min()*1e3),float(h.max()*1e3),float(h.max()/h.min()),loc,demand,(time.perf_counter()-t0)*1000)
                rows.append(row)
                if n==49 and strategy==strategies[-1]:
                    eta,nov,locality,E2,weights=pulse_data(model,U,scenario.freqs)
                    maps[scenario.name]={'nodes_m':x.tolist(),'indicator':normalize(eta).tolist(),'novelty':nov.tolist(),'locality_index':locality,'coverage_demand':demand}
    elapsed=time.perf_counter()-start
    return rows,maps,elapsed


def summarize(rows, cycle):
    summary={'cycle':cycle,'by_scenario':{},'aggregate':{}}
    scenarios=sorted(set(r.scenario for r in rows)); strategies=sorted(set(r.strategy for r in rows))
    for sc in scenarios:
        summary['by_scenario'][sc]={}
        for st in strategies:
            rs=[r for r in rows if r.scenario==sc and r.strategy==st]
            summary['by_scenario'][sc][st]={
                'mean_qoi_error':float(np.mean([r.qoi_rel_error for r in rs])),
                'mean_field_error':float(np.mean([r.field_rel_error for r in rs])),
                'qoi_error_49':next((r.qoi_rel_error for r in rs if r.nodes==49),None),
                'field_error_49':next((r.field_rel_error for r in rs if r.nodes==49),None),
                'best_qoi_error':min(r.qoi_rel_error for r in rs),
                'mean_h_ratio':float(np.mean([r.h_ratio for r in rs])),
            }
    for st in strategies:
        rs=[r for r in rows if r.strategy==st]
        summary['aggregate'][st]={
            'geometric_mean_qoi_error':float(np.exp(np.mean(np.log([r.qoi_rel_error+1e-15 for r in rs])))),
            'geometric_mean_field_error':float(np.exp(np.mean(np.log([r.field_rel_error+1e-15 for r in rs])))),
            'mean_runtime_ms':float(np.mean([r.build_runtime_ms for r in rs])),
            'mean_h_ratio':float(np.mean([r.h_ratio for r in rs])),
        }
    return summary


def write_cycle(cycle,rows,maps,elapsed,summary):
    stem=OUT/f'piezo_evolution_cycle{cycle}'
    with open(stem.with_suffix('.csv'),'w',newline='') as fp:
        w=csv.DictWriter(fp,fieldnames=list(asdict(rows[0]).keys())); w.writeheader(); w.writerows(asdict(r) for r in rows)
    payload={'elapsed_s':elapsed,'summary':summary,'maps':maps,'rows':[asdict(r) for r in rows]}
    with open(stem.with_suffix('.json'),'w') as fp: json.dump(payload,fp,indent=2)
    with open(stem.with_name(stem.name+'_summary.txt'),'w') as fp:
        fp.write(f'Cycle {cycle} elapsed_s={elapsed:.3f}\n')
        for sc,data in summary['by_scenario'].items():
            fp.write(f'\n{sc}\n')
            for st,v in data.items():
                fp.write(f"  {st:20s} qoi49={100*v['qoi_error_49']:.4f}% field49={100*v['field_error_49']:.4f}% mean_qoi={100*v['mean_qoi_error']:.4f}% hratio={v['mean_h_ratio']:.2f}\n")
        fp.write('\nAggregate geometric means\n')
        for st,v in summary['aggregate'].items():
            fp.write(f"  {st:20s} qoi={100*v['geometric_mean_qoi_error']:.4f}% field={100*v['geometric_mean_field_error']:.4f}% runtime={v['mean_runtime_ms']:.2f}ms hratio={v['mean_h_ratio']:.2f}\n")


def main():
    scenarios=make_scenarios()
    cycle1_strategies=('uniform','snapshot','pulse_v0','cycle1_regularized')
    rows1,maps1,t1=run_cycle(1,cycle1_strategies,scenarios)
    sum1=summarize(rows1,1); write_cycle(1,rows1,maps1,t1,sum1)

    cycle2_strategies=('uniform','snapshot','pulse_v0','cycle1_regularized','cycle2_governed')
    rows2,maps2,t2=run_cycle(2,cycle2_strategies,scenarios)
    sum2=summarize(rows2,2); write_cycle(2,rows2,maps2,t2,sum2)

    combined={
        'method':'Pulsed Residual-Governed Sparsening — two-cycle evolution',
        'cycle1':{'elapsed_s':t1,'summary':sum1},
        'cycle2':{'elapsed_s':t2,'summary':sum2},
        'scenario_descriptions':{s.name:s.description for s in scenarios},
        'mathematical_notes':[
            'Linear quasi-static scalar load progression remains rank-identical and is not used as evidence of dynamic adaptivity.',
            'Cycle 1 adds wave-resolution pressure, mesh-size coverage, neighbor smoothing, and bulk marking.',
            'Cycle 2 computes spatial entropy of pulse-integrated residual energy and uses it as an authority coefficient between a localized residual channel and a global coverage channel.',
            'The compression gate closes continuously as effective residual support approaches the full domain.'
        ]
    }
    with open(OUT/'piezo_evolution_two_cycle_results.json','w') as fp: json.dump(combined,fp,indent=2)
    print((OUT/'piezo_evolution_cycle1_summary.txt').read_text())
    print((OUT/'piezo_evolution_cycle2_summary.txt').read_text())


if __name__=='__main__':
    main()
