"""Generates sim/fig4_data.pkl: GLS trajectory sweep (20 realizations per Na), envelope ratios,
full RMS(z) curves, one representative Na=1e7 realization, and the clean-pipeline reference.
Requires: sim/stage4_design.pkl (from stage4_gls.py), sim/true_traj.pkl (from stage2_run.py).
Runtime ~2 min."""
import numpy as np, pickle, sys, time
sys.path.insert(0,'sim')
import stage4_gls as G, stage5_alloc as S
TT=pickle.load(open('sim/true_traj.pkl','rb'))
zeval=np.linspace(0.05,0.60,60)
out={'zeval':zeval,'TT_final_std':TT[:,-1].std()}
for Na in [1e4,1e5,1e6,1e7]:
    t0=time.time(); curves=[]; envs=[]; keep=None
    for r in range(20):
        rng=np.random.default_rng(70000+r)
        rows=[]
        for k in range(20):
            P=S.DES[k]
            mom,sg=G.mom_noisy(P,rng,int(Na))
            rows.append(G.gls_fit(P,mom-P['rw'],sg))
        V=S.resample(rows); T=S.traj(V)
        curves.append(np.sqrt(np.mean((T-TT)**2,axis=0)))
        envs.append(T[:,-1].std()/TT[:,-1].std())
        if r==0 and Na==1e7: keep=(V.astype(np.float32),T.astype(np.float32))
    curves=np.array(curves); envs=np.array(envs)
    out[Na]=dict(rms_curves=curves.astype(np.float32),env=envs)
    if keep: out['rep_1e7']=keep
    print(f"Na={Na:.0e}: {curves[:,-1].mean()*1e6:.1f} um  [{time.time()-t0:.0f}s]")
rows=[G.gls_fit(S.DES[k],S.DES[k]['mom_exact']-S.DES[k]['rw'],np.ones(S.DES[k]['xc'].size)) for k in range(20)]
Vc=S.resample(rows); Tc=S.traj(Vc)
out['clean']=dict(V=Vc.astype(np.float32),rms=np.sqrt(np.mean((Tc-TT)**2,axis=0)),env=Tc[:,-1].std()/TT[:,-1].std(),T=Tc.astype(np.float32))
pickle.dump(out,open('sim/fig4_data.pkl','wb')); print("saved sim/fig4_data.pkl")
