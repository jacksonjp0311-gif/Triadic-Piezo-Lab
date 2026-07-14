from __future__ import annotations
import csv,json,time,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,'/mnt/data')
import piezo_evolution_two_cycle as core
import piezo_triad_deep_benchmark as tri
import piezo_triadic_stability_bridge as bridge

OUT=Path('/mnt/data')


def certificate(metrics, transfer_tol=.003, phase_tol_deg=.25, point_tol=.02):
    channels={
      'transfer':metrics['transfer_rel_error']/transfer_tol,
      'phase':metrics['phase_rms_deg']/phase_tol_deg,
      'point':metrics['max_regularized_point_error']/point_tol,
    }
    return {'channels':channels,'delta_stability':float(max(channels.values())),
            'pass':bool(max(channels.values())<=1.0)}


def enrich_multigate(ctx,freqs,Ufine,transfer_tol=.003,phase_tol_deg=.25,point_tol=.02,max_modes=72,batch=2):
    ncomp=len(ctx.components);q=np.zeros(ncomp,int);active=bridge.active_indices(len(freqs),5);history=[]
    adj=np.zeros((len(freqs),len(ctx.model.x)-1),complex);ctip=np.zeros(len(ctx.model.x)-1,complex);ctip[-1]=1
    for j,f in enumerate(freqs):
        lo,di,up=tri.system_diagonals(ctx.model,f)
        adj[j]=tri.thomas(np.conjugate(up),np.conjugate(di),np.conjugate(lo),ctip)
    U,meta=bridge.solve_cb(ctx,q,freqs);met=bridge.transfer_metrics(U,Ufine,freqs);cert=certificate(met,transfer_tol,phase_tol_deg,point_tol)
    while not cert['pass'] and int(q.sum())<max_modes:
        hf=Ufine[:,-1];hr=U[:,-1]
        point=np.abs(hr-hf)/(np.abs(hf)+.02*np.max(np.abs(hf))+1e-30)
        phase=np.abs(np.angle(hr*np.conjugate(hf)))
        # Promote the frequency that most violates either point or phase coherence.
        violation=np.maximum(point/point_tol,phase/(phase_tol_deg*np.pi/180+1e-30))
        worst=int(np.argmax(violation))
        if worst not in active:active.append(worst);active.sort()
        scores=np.full(ncomp,-np.inf,float)
        for c,comp in enumerate(ctx.components):
            k=int(q[c])
            if k>=comp.phi.shape[1] or len(comp.interior_global)==0:continue
            v=comp.phi[:,k];rows=comp.interior_global-1;score=0.0
            for j in active:
                f=freqs[j];w=2*np.pi*f;Z=ctx.full_K-w*w*ctx.full_M+1j*w*ctx.full_C
                r=ctx.full_F-Z@U[j,1:]
                rv=np.vdot(v,r[rows]);zv=np.vdot(v,adj[j,rows])
                den=abs(np.vdot(v,Z[np.ix_(rows,rows)]@v))+1e-30
                ampw=.25+.75*abs(Ufine[j,-1])/(np.max(np.abs(Ufine[:,-1]))+1e-30)
                # Frequencies currently violating phase/point gates receive authority.
                vw=1+min(5.0,float(violation[j]))
                score+=ampw*vw*abs(rv)*abs(zv)/den
            scores[c]=score/(1+.02*k)
        available=np.where(np.isfinite(scores))[0]
        if not len(available):break
        take=min(batch,max_modes-int(q.sum()),len(available))
        chosen=available[np.argsort(scores[available])[::-1][:take]]
        if np.max(scores[available])<=1e-30:
            cand=[(ctx.components[c].eigfreqs[q[c]],c) for c in available]
            chosen=np.array([c for _,c in sorted(cand)[:take]],int)
        for c in chosen:q[c]+=1
        U,meta=bridge.solve_cb(ctx,q,freqs);met=bridge.transfer_metrics(U,Ufine,freqs);cert=certificate(met,transfer_tol,phase_tol_deg,point_tol)
        history.append({'step':len(history)+1,'chosen_components':[int(c) for c in chosen],
                        'q_alloc':q.tolist(),'active_frequency_count':len(active),
                        'active_frequencies_hz':[float(freqs[i]) for i in active],**met,**cert})
    return q,U,meta,met,cert,history,active


def eval_one(sc,fine_nodes=193,stride=16,ref_nodes=513,transfer_tol=.003,phase_tol=.25,point_tol=.02,max_modes=72):
    x=np.linspace(0,core.L,fine_nodes);m=core.assemble(x,sc);Ufine=tri.solve_harmonic_fast(m,sc.freqs);ctx=bridge.build_cb_context(m,stride)
    t=time.perf_counter();q,U,meta,met,cert,hist,active=enrich_multigate(ctx,sc.freqs,Ufine,transfer_tol,phase_tol,point_tol,max_modes)
    elapsed=(time.perf_counter()-t)*1000
    qspec=bridge.spectral_allocation(ctx,int(q.sum()));Us,ms=bridge.solve_cb(ctx,qspec,sc.freqs);mets=bridge.transfer_metrics(Us,Ufine,sc.freqs);certs=certificate(mets,transfer_tol,phase_tol,point_tol)
    stability=bridge.stability_certificate(meta)
    return {'scenario':sc.name,'fmin':float(sc.freqs.min()),'fmax':float(sc.freqs.max()),'full_dofs':fine_nodes-1,'port_dofs':len(ctx.port_col),
      'multigate':{'q_alloc':q.tolist(),'retained_modes':int(q.sum()),'reduced_dofs':meta['reduced_dofs'],'reduction':meta['reduction'],
                   'runtime_ms':elapsed,'active_count':len(active),'history':hist,**met,**cert,'stability':stability},
      'spectral_same_modes':{'q_alloc':qspec.tolist(),'retained_modes':int(qspec.sum()),'reduced_dofs':ms['reduced_dofs'],'reduction':ms['reduction'],**mets,**certs,
                             'stability':bridge.stability_certificate(ms)}}


def gm(x):return float(np.exp(np.mean(np.log(np.maximum(np.asarray(x,float),1e-18)))))

def main():
    det=[]
    for sc in core.make_scenarios():
        r=eval_one(sc);det.append(r);m=r['multigate'];print('DET',sc.name,m['pass'],m['retained_modes'],m['transfer_rel_error'],m['phase_rms_deg'],m['max_regularized_point_error'])
    rng=np.random.default_rng(20260715);rnd=[]
    for i in range(24):
        sc,meta=tri.random_scenario(rng,i);r=eval_one(sc);r['random_meta']=meta;rnd.append(r);m=r['multigate'];print('RND',i,m['pass'],m['retained_modes'],m['transfer_rel_error'],m['phase_rms_deg'],m['max_regularized_point_error'])
    rows=det+rnd
    mm=[r['multigate'] for r in rows];ss=[r['spectral_same_modes'] for r in rows]
    summary={'count':len(rows),'passes':sum(m['pass'] for m in mm),'median_modes':float(np.median([m['retained_modes'] for m in mm])),
      'mean_modes':float(np.mean([m['retained_modes'] for m in mm])),'median_reduction':float(np.median([m['reduction'] for m in mm])),
      'geomean_transfer_error':gm([m['transfer_rel_error'] for m in mm]),'median_phase_deg':float(np.median([m['phase_rms_deg'] for m in mm])),
      'max_phase_deg':float(np.max([m['phase_rms_deg'] for m in mm])),'max_point_error':float(np.max([m['max_regularized_point_error'] for m in mm])),
      'vs_spectral_wins_delta':sum(m['delta_stability']<s['delta_stability'] for m,s in zip(mm,ss)),
      'spectral_passes_same_modes':sum(s['pass'] for s in ss),'all_stability_positive':all(m['stability']['positive_K'] and m['stability']['positive_M'] and m['stability']['nonnegative_C'] for m in mm),
      'law':'Compression is stable only when global transfer, accumulated phase, and extremal frequency mismatch pass simultaneously on an invariant port manifold.'}
    payload={'certificate_tolerances':{'transfer':.003,'phase_deg':.25,'point':.02},'summary':summary,'deterministic':det,'random_holdout':rnd}
    (OUT/'piezo_triadic_multigate_results.json').write_text(json.dumps(payload,indent=2))
    flat=[]
    for kind,rs in [('deterministic',det),('random',rnd)]:
      for r in rs:
       m=r['multigate'];s=r['spectral_same_modes'];flat.append({'kind':kind,'scenario':r['scenario'],'fmax':r['fmax'],'modes':m['retained_modes'],'dofs':m['reduced_dofs'],'reduction':m['reduction'],'pass':m['pass'],'delta_stability':m['delta_stability'],'transfer_error':m['transfer_rel_error'],'phase_deg':m['phase_rms_deg'],'point_error':m['max_regularized_point_error'],'spectral_pass':s['pass'],'spectral_delta':s['delta_stability']})
    with open(OUT/'piezo_triadic_multigate.csv','w',newline='') as fp:
      w=csv.DictWriter(fp,fieldnames=flat[0].keys());w.writeheader();w.writerows(flat)
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
