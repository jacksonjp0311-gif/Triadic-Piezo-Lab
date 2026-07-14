from __future__ import annotations
import csv,json,sys,time
from pathlib import Path
import numpy as np
sys.path.insert(0,'/mnt/data')
import piezo_axisymmetric_2d_evolution as a
import piezo_axisymmetric_2d_port_evolution as pe
from piezo_axisymmetric_2d_adversarial_closure import point_violation
OUT=Path('/mnt/data')

def evaluate(sc,level=2,p_target=16,maxadd=8):
    t=time.perf_counter();m=a.assemble(sc,level);Uf=a.solve_full(m,sc.freqs);base=a.build_components(m,max_modes_per=16)
    train=sorted(set(int(round(x)) for x in np.linspace(0,len(sc.freqs)-1,7)));hist=[]
    for it in range(maxadd+1):
        pc,en,s=pe.port_context(base,Uf[train],p_target)
        rr=a.enrich(pc,sc.freqs[train],Uf[train],2,max_modes=180,batch=4)
        U,meta=a.solve_rom(pc,rr['q'],sc.freqs);met=a.all_metrics(m,U,Uf);delta,ch=a.cert_cycle2(met)
        hist.append({'iteration':it,'sentinel_indices':train.copy(),'sentinel_count':len(train),'port_modes':len(pc.ports),'local_modes':int(rr['q'].sum()),'global_dofs':meta['dofs'],'reduction':meta['reduction'],'delta':float(delta),'channels':ch,'metrics':met,'training_delta':rr['delta']})
        if delta<=1:break
        vv=point_violation(m,U,Uf);vv[train]=-np.inf;w=int(np.argmax(vv))
        if not np.isfinite(vv[w]):break
        train.append(w);train.sort()
    last=hist[-1]
    return {'scenario':sc.name,'description':sc.description,'level':level,'full_dofs':m.K.shape[0],'elements':len(m.elements),'nodal_ports':len(base.ports),'pass':bool(last['delta']<=1),'sentinels':last['sentinel_count'],'sentinel_indices':last['sentinel_indices'],'sentinel_frequencies':[float(sc.freqs[i]) for i in last['sentinel_indices']],'port_modes':last['port_modes'],'local_modes':last['local_modes'],'global_dofs':last['global_dofs'],'reduction':last['reduction'],'delta':last['delta'],'channels':last['channels'],'metrics':last['metrics'],'history':hist,'elapsed_s':time.perf_counter()-t}

def gm(x):return float(np.exp(np.mean(np.log(np.maximum(np.asarray(x,float),1e-18)))))
def main():
    rows=[]
    for sc in a.scenario_set():
        r=evaluate(sc);r['kind']='deterministic';rows.append(r);print('DET',sc.name,r['pass'],'sent',r['sentinels'],'p',r['port_modes'],'m',r['local_modes'],'dof',r['global_dofs'],'red',round(r['reduction'],4),'delta',round(r['delta'],3),flush=True)
    rng=np.random.default_rng(20260714)
    for i in range(10):
        r=evaluate(a.random_scenario(rng,i));r['kind']='holdout';rows.append(r);print('RND',i,r['pass'],'sent',r['sentinels'],'p',r['port_modes'],'m',r['local_modes'],'dof',r['global_dofs'],'red',round(r['reduction'],4),'delta',round(r['delta'],3),flush=True)
    summary={'count':len(rows),'passes':sum(r['pass'] for r in rows),'median_full_dofs':float(np.median([r['full_dofs'] for r in rows])),'median_nodal_ports':float(np.median([r['nodal_ports'] for r in rows])),'median_sentinels':float(np.median([r['sentinels'] for r in rows])),'max_sentinels':int(max(r['sentinels'] for r in rows)),'median_port_modes':float(np.median([r['port_modes'] for r in rows])),'median_local_modes':float(np.median([r['local_modes'] for r in rows])),'median_global_dofs':float(np.median([r['global_dofs'] for r in rows])),'median_reduction':float(np.median([r['reduction'] for r in rows])),'min_reduction':float(min(r['reduction'] for r in rows)),'max_delta':float(max(r['delta'] for r in rows)),'geomean_transfer':gm([r['metrics']['transfer'] for r in rows]),'geomean_surface':gm([r['metrics']['surface_rel'] for r in rows]),'geomean_energy':gm([r['metrics']['energy_rel'] for r in rows]),'law':'As the physical manifold is refined, globally exposed complexity need not grow proportionally: a frequency-closed port basis plus nested local modes can preserve the fine 2D continuum while the global system remains small.'}
    payload={'method':'Fine 2D axisymmetric adversarially closed triadic port bridge','summary':summary,'rows':rows}
    (OUT/'piezo_axisymmetric_2d_fine_bridge_results.json').write_text(json.dumps(payload,indent=2))
    flat=[]
    for r in rows:
        m=r['metrics'];flat.append({'kind':r['kind'],'scenario':r['scenario'],'pass':r['pass'],'full_dofs':r['full_dofs'],'elements':r['elements'],'nodal_ports':r['nodal_ports'],'sentinels':r['sentinels'],'port_modes':r['port_modes'],'local_modes':r['local_modes'],'global_dofs':r['global_dofs'],'reduction':r['reduction'],'delta':r['delta'],'transfer':m['transfer'],'phase_deg':m['phase_deg'],'point':m['point'],'surface':m['surface_rel'],'shape':m['shape'],'energy':m['energy_rel'],'energy_distribution':m['energy_distribution']})
    with open(OUT/'piezo_axisymmetric_2d_fine_bridge_results.csv','w',newline='') as fp:
        w=csv.DictWriter(fp,fieldnames=flat[0]);w.writeheader();w.writerows(flat)
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
