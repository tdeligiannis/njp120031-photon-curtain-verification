"""Stage 2 reference -- GATE ONLY (v8 restructure, de-hazarded).
v7 note: this script formerly built its own 251-bin maps and 30 realizations, colliding
with stage2_maps.pkl / stage2_results.pkl formats. In v8 it loads the canonical 2500-bin
adaptive maps (from stage2_maps.py) and runs Gate 1: the cumulative-estimator identity
(Eq. 24 path) against the direct formula (Eq. 13 path), plane by plane."""
import numpy as np, pickle, sys
sys.path.insert(0,'sim')
import v8_config as C
from sim_core import cum_integrate
W=C.W
MAPS=pickle.load(open('sim/stage2_maps.pkl','rb'))
NB=MAPS[0]['rho_b'].size; xf_edges=np.linspace(-1.25e-3,1.25e-3,NB+1); xf=0.5*(xf_edges[1:]+xf_edges[:-1])
print("GATE 1 (v8): cumulative-vs-direct v, max rel dev on valid & |v_f|>2 um/s")
out={}
for k,mp in enumerate(MAPS):
    prodb=(mp['a_w']*mp['rho_b'][None,:]).astype(float)
    lev=xf[None,:]-mp['xc'][:,None]
    mom=((lev*prodb).sum(axis=1)+mp['D_corr'])/W
    vcum=cum_integrate(mom-mp['rho_w'],mp['xc'])/(mp['t']*np.maximum(mp['rho_w'],1e-300))
    sel=mp['valid']&(np.abs(mp['v_formula'])>2e-6)
    dev=np.max(np.abs((vcum[sel]-mp['v_formula'][sel])/mp['v_formula'][sel]))
    out[k]=dict(z=mp['z'],dev=dev,npts=int(sel.sum()))
    if k in [0,1,9,19] or dev>0.05:
        print(f"  z={mp['z']*100:5.1f} cm : {dev*100:6.2f}%   ({sel.sum()} pts)")
devs=np.array([out[k]['dev'] for k in range(len(MAPS))])
print(f"GATE 1 worst plane: {devs.max()*100:.2f}%   median plane: {np.median(devs)*100:.2f}%")
pickle.dump(out,open('sim/stage2ref_gate.pkl','wb'))
print("PASS" if devs.max()<0.05 else "FAIL", "- saved sim/stage2ref_gate.pkl")
