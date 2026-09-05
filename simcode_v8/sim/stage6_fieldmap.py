"""Quantity 1: early-plane field-map campaign. Planes {0,1,2} (z=5.0/7.7/10.5cm).
Metrics per (plane, Na): core-region relative & absolute field error (12 realizations),
matched-filter SNR (total + oscillatory fringe-channeling component)."""
import numpy as np, pickle, os, sys, time
sys.path.insert(0,'sim')
import stage4_gls as G
from scipy.ndimage import gaussian_filter1d
import v8_config as C
W=C.W
zj=np.linspace(0.05,0.57,20)
DES=pickle.load(open('sim/stage4_design.pkl','rb'))
Vt=[G.gls_fit(DES[k],DES[k]['mom_exact']-DES[k]['rw'],np.ones(DES[k]['xc'].size)) for k in range(20)]
R4=pickle.load(open('sim/stage4_results.pkl','rb')); snr0=R4['snr_1e4']
def sig_of(P,Na):
    lev=G.xf[None,:]-P['xc'][:,None]
    return np.sqrt(((lev**2)*P['rho']*P['win']).sum(axis=1)*P['SIG']**2/(Na*W**2))
PLANES=[0,1,2]; NAS=[3e5,1e6,3e6,1e7,3e7]; NREAL=12
fn='sim/stage6_results.pkl'
store=pickle.load(open(fn,'rb')) if os.path.exists(fn) else {}
# deterministic: oscillatory-component SNR per plane
if 'osc' not in store:
    store['osc']={}
    for k in PLANES:
        P=DES[k]; xc=P['xc']; lam_f=C.lam_f(zj[k])
        sm=gaussian_filter1d(Vt[k],sigma=(lam_f/2)/(xc[1]-xc[0]),mode='nearest')
        vosc=Vt[k]-sm
        mom_osc=P['t']*np.gradient(P['rw']*vosc,xc)
        sg4=sig_of(P,1e4)
        snr_osc=np.sqrt(((mom_osc/sg4)**2).sum())
        store['osc'][k]=dict(snr_osc_1e4=snr_osc,vscale=np.median(np.abs(Vt[k][P['rw']>=0.05*P['rw'].max()])),
                             osc_amp=np.sqrt(np.mean(vosc[P['rw']>=0.05*P['rw'].max()]**2)))
        print(f"plane {k} (z={zj[k]*100:.1f}cm): field scale {store['osc'][k]['vscale']*1e6:.1f} um/s, "
              f"osc amp {store['osc'][k]['osc_amp']*1e6:.1f} um/s, SNR_osc(1e4)={snr_osc:.3f}, SNR_tot(1e4)={snr0[k]:.3f}")
    pickle.dump(store,open(fn,'wb'))
for k in PLANES:
    P=DES[k]; core=P['rw']>=0.05*P['rw'].max()
    vs=np.median(np.abs(Vt[k][core]))
    for Na in NAS:
        key=(k,Na)
        if key in store: continue
        t0=time.time(); rel=[]; ab=[]
        for r in range(NREAL):
            rng=np.random.default_rng(50000+1000*k+r)
            mom,sg=G.mom_noisy(P,rng,int(Na))
            vg=G.gls_fit(P,mom-P['rw'],sg)
            e=np.sqrt(np.mean((vg-Vt[k])[core]**2))
            ab.append(e); rel.append(e/vs)
        store[key]=dict(rel=np.array(rel),abs=np.array(ab))
        pickle.dump(store,open(fn,'wb'))
        print(f"  plane {k}, Na={Na:.0e}: rel err {np.mean(rel)*100:6.2f} +- {np.std(rel)/np.sqrt(NREAL)*100:4.2f}%   "
              f"abs {np.mean(ab)*1e6:7.2f} um/s   [{time.time()-t0:.0f}s]")
print("\nDONE — sim/stage6_results.pkl")
