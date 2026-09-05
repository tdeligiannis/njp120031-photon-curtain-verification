"""Stage 5 / Quantity 2: atom-allocation optimization across planes.
Per-plane trajectory-error coefficients A_k by single-plane noise injection through GLS pipeline.
sigma_traj^2 = sum_k A_k/N_k (gated by additivity + 1/Na scaling); optimal N_k ~ sqrt(A_k).
Usage: python3 sim/stage5_alloc.py <Na> <nreal>"""
import numpy as np, pickle, os, sys, time
sys.path.insert(0,'sim')
import stage4_gls as G
from scipy.interpolate import RectBivariateSpline
from sim_core import make_psi
import v8_config as C
d,sig0,L=C.d,C.sig0,C.L; W=C.W
zj=np.linspace(0.05,0.57,20)
DES=pickle.load(open('sim/stage4_design.pkl','rb'))
psi,_,_=make_psi(d,sig0,tmap=C.t_of_z)
xr=np.linspace(-950e-6,950e-6,381)
xgf=np.linspace(-3e-3,3e-3,2**15); p0=np.abs(psi(xgf,0.05))**2; cdf=np.cumsum(p0)/p0.sum()
x0=np.interp((np.arange(200)+0.5)/200,cdf,xgf); zeval=np.linspace(0.05,L,60)
def traj(V2):
    sp=RectBivariateSpline(zj,xr,V2,kx=3,ky=3)
    f=lambda z,x: sp(z,np.clip(x,xr[0],xr[-1]),grid=False)/C.vz_of(z)
    zf=np.linspace(0.05,L,221); h=zf[1]-zf[0]; x=x0.copy(); X=[x.copy()]
    for i in range(zf.size-1):
        k1=f(zf[i],x);k2=f(zf[i]+h/2,x+h/2*k1);k3=f(zf[i]+h/2,x+h/2*k2);k4=f(zf[i]+h,x+h*k3)
        x=x+h/6*(k1+2*k2+2*k3+k4); X.append(x.copy())
    return np.array([np.interp(zeval,zf,np.array(X).T[j]) for j in range(200)])
def resample(rows):
    out=np.zeros((20,xr.size))
    for k in range(20):
        P=DES[k]; out[k]=np.interp(xr,P['xc'],rows[k],left=0,right=0)
        out[k][(xr<P['lo'])|(xr>P['hi'])]=0
    return out
Vt=[G.gls_fit(DES[k],DES[k]['mom_exact']-DES[k]['rw'],np.ones(DES[k]['xc'].size)) for k in range(20)]
Tref=traj(resample(Vt))

def run_alloc(shares,Ntot,nreal,seed0):
    rmss=[]
    for r in range(nreal):
        rng=np.random.default_rng(seed0+r)
        rows=[]
        for j in range(20):
            P=DES[j]; Nj=max(int(Ntot*shares[j]),200)
            mom,sg=G.mom_noisy(P,rng,Nj)
            rows.append(G.gls_fit(P,mom-P['rw'],sg))
        T=traj(resample(rows))
        rmss.append(np.sqrt(np.mean((T[:,-1]-Tref[:,-1])**2)))
    return np.array(rmss)

if __name__=="__main__":
    Na=int(float(sys.argv[1])); NREAL=int(sys.argv[2])
    fn=f'sim/stage5_inject_Na{Na:.0e}.pkl'.replace('+0','')
    store=pickle.load(open(fn,'rb')) if os.path.exists(fn) else {}
    # zeroth-order estimate from stored matched-filter SNRs (uniform weights)
    if 'zeroth' not in store:
        R4=pickle.load(open('sim/stage4_results.pkl','rb')); snr=R4['snr_1e4']
        c=1/snr**2
        store['zeroth']=np.sqrt(20*c.sum())/np.sqrt(c).sum()
        print(f"zeroth-order gain (uniform w, A_k ~ 1/SNR_k^2): G = {store['zeroth']:.2f}")
    t0=time.time()
    def sig_of(P,Na):
        lev=G.xf[None,:]-P['xc'][:,None]
        var=((lev**2)*P['rho']*P['win']).sum(axis=1)*P['SIG']**2/(Na*W**2)
        return np.sqrt(var)
    for k in range(20):
        key=f'A_{k}'
        if key in store: continue
        P=DES[k]; sgk=sig_of(P,Na)
        # bias-matched reference: same weights+ridge on exact data
        ref_k=G.gls_fit(P,P['mom_exact']-P['rw'],sgk)
        rows_ref=[Vt[j] for j in range(20)]; rows_ref[k]=ref_k
        Tref_k=traj(resample(rows_ref))
        rr=[]
        for r in range(NREAL):
            rng=np.random.default_rng(9000+100*k+r)
            mom,sg=G.mom_noisy(P,rng,Na)
            rows=[Vt[j] for j in range(20)]; rows[k]=G.gls_fit(P,mom-P['rw'],sg)
            T=traj(resample(rows))
            rr.append(np.mean((T[:,-1]-Tref_k[:,-1])**2))
        store[key]=dict(msq=np.mean(rr),se=np.std(rr)/np.sqrt(NREAL))
        pickle.dump(store,open(fn,'wb'))
        print(f"  plane {k:2d} (z={zj[k]*100:5.1f}cm): RMS_k = {np.sqrt(store[key]['msq'])*1e6:8.1f} um   [{time.time()-t0:.0f}s]")
    A=np.array([store[f'A_{k}']['msq'] for k in range(20)])*Na
    store['A']=A; pickle.dump(store,open(fn,'wb'))
    print(f"\nadditivity check: sqrt(sum A_k/Na) = {np.sqrt(A.sum()/Na)*1e6:.1f} um  (compare full-noise GLS RMS at this Na)")
    G_gain=np.sqrt(20*A.sum())/np.sqrt(A).sum()
    print(f"optimal-allocation gain (this Na, additive model): G = {G_gain:.3f}")
    print("optimal shares N_k/N_tot:", " ".join(f"{s:.3f}" for s in np.sqrt(A)/np.sqrt(A).sum()))
    # ---- m16: seed-paired (common-random-number) allocation comparison ----
    shares_opt=np.sqrt(A)/np.sqrt(A).sum(); shares_uni=np.ones(20)/20.0
    NPAIR=12; s0=777000
    r_uni=run_alloc(shares_uni,20*Na,NPAIR,s0)
    r_opt=run_alloc(shares_opt,20*Na,NPAIR,s0)   # identical seed list: CRN pairing
    gain=r_uni/r_opt
    store['m16']=dict(Na=Na,shares_opt=shares_opt,r_uni=r_uni,r_opt=r_opt,gain=gain)
    pickle.dump(store,open(fn,'wb'))
    print(f"m16 CRN-paired allocation: uniform/optimal RMS ratio = {gain.mean():.3f} +- "
          f"{gain.std(ddof=1)/np.sqrt(NPAIR):.3f}   ({(gain.mean()-1)*100:+.1f}% benefit)")

