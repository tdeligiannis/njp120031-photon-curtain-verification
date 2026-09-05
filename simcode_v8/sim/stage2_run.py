"""Stage 2 runner: 30 noise realizations at reference cell (eta*eps=0), Na=1e4.
Restartable: python3 sim/stage2_run.py <n_target>. Checkpoints every 5 realizations."""
import numpy as np, pickle, os, sys, time
sys.path.insert(0,'sim'); from sim_core import make_psi, hbar, M
import v8_config as C
from scipy.interpolate import RectBivariateSpline
from scipy.integrate import solve_ivp

d,sig0,L = C.d,C.sig0,C.L; W=C.W
Na=int(1e4)   # v8: per-plane SIG/V/Psc carried in the maps; SE bg = Psc(z) x recoil-smeared rho
MAPS=pickle.load(open('sim/stage2_maps.pkl','rb'))
NB=MAPS[0]['rho_b'].size; xf_edges=np.linspace(-1.25e-3,1.25e-3,NB+1); xf=0.5*(xf_edges[1:]+xf_edges[:-1])
zj=np.array([mp['z'] for mp in MAPS])
psi,dpsi,v_true=make_psi(d,sig0,tmap=C.t_of_z)
xr=np.linspace(-950e-6,950e-6,381)

from sim_core import cum_integrate

def vfield(mp,rng=None):
    """One plane's reconstructed v on xr. Windowed data moment + plug-in outside/sub-bin terms.
    Each x_c is an independent ensemble: counts drawn per curtain position."""
    rhoV=mp['V']*mp['rho_b']+(1-mp['V'])*mp['rho_b_inc']            # m12 visibility mixture
    base=(1-mp['Psc'])*rhoV+mp['Psc']*mp['bg_b']                    # per-plane SE background
    sh=lambda s: np.interp(xf,xf+s,base,left=0,right=0)
    pmix=(1-C.F_IMP)*base+0.5*C.F_IMP*(sh(+C.DSCR_IMP)+sh(-C.DSCR_IMP))   # f-impurity split
    rho=mp['rho_b']                                  # ideal coherent density: analyst-side model
    prodb=mp['a_w']*rho[None,:]                      # model rho*a per (xc,bin)
    recw=(1-C.F_IMP)*(1-mp['Psc'])*(mp['V']*prodb
         +(1-mp['V'])*mp['a_w_inc'].astype(float)*mp['rho_b_inc'][None,:])   # mixture mean record-weighted density
    ahm=np.where(pmix[None,:]>0,recw/np.maximum(pmix[None,:],1e-300),0.0)    # mean record per detected atom
    win=np.abs(prodb)>1e-3*np.abs(prodb).max(axis=1,keepdims=True)
    lev=xf[None,:]-mp['xc'][:,None]
    Mout=(lev*prodb*(~win)).sum(axis=1)              # deterministic outside-window moment
    if rng is None:
        dat=(lev*recw*win).sum(axis=1)
        dat/= (1-C.F_IMP)*(1-mp['Psc'])   # analyst normalization: F and Psc independently known (purity spec, photon budget)
    else:
        pvec=np.concatenate([pmix,[max(1-pmix.sum(),0)]]); pvec/=pvec.sum()
        n=rng.multinomial(Na,pvec,size=mp['xc'].size)[:,:NB]      # per-x_c ensembles
        ah=ahm+rng.standard_normal((mp['xc'].size,NB))*mp['SIG']/np.sqrt(np.maximum(n,1))
        dat=(lev*(n/Na)*ah*win).sum(axis=1)
        dat/= (1-C.F_IMP)*(1-mp['Psc'])   # analyst normalization (same known factors)
    mom=(dat+Mout+mp['D_corr'])/W
    rw_use=mp['V']*mp['rho_w']+(1-mp['V'])*mp['rho_w_inc']   # self-consistent measured window density (Sec.-7 calibration semantics)
    numc=cum_integrate(mom-rw_use,mp['xc'])
    v=numc/(mp['t']*np.maximum(rw_use,1e-300)); v[~mp['valid']]=0.0
    vr=np.interp(xr,mp['xc'],v,left=0.0,right=0.0)
    vr[(xr<mp['xc'][mp['valid']].min())|(xr>mp['xc'][mp['valid']].max())]=0.0
    return vr

NTRAJ=200; zs0=0.05
xgf=np.linspace(-3e-3,3e-3,2**15); p0=np.abs(psi(xgf,zs0))**2
cdf=np.cumsum(p0)/p0.sum(); x0=np.interp((np.arange(NTRAJ)+0.5)/NTRAJ,cdf,xgf)
zeval=np.linspace(zs0,L,60)
def traj_recon(V):
    sp=RectBivariateSpline(zj,xr,V,kx=3,ky=3)
    f=lambda z,x: sp(z,np.clip(x,xr[0],xr[-1]),grid=False)/C.vz_of(z)
    zfine=np.linspace(zs0,L,221); h=zfine[1]-zfine[0]
    X=np.empty((x0.size,zfine.size)); X[:,0]=x0; x=x0.copy()
    for i in range(zfine.size-1):
        z=zfine[i]
        k1=f(z,x); k2=f(z+h/2,x+h/2*k1); k3=f(z+h/2,x+h/2*k2); k4=f(z+h,x+h*k3)
        x=x+h/6*(k1+2*k2+2*k3+k4); X[:,i+1]=x
    return np.array([np.interp(zeval,zfine,X[j]) for j in range(x0.size)])
if not os.path.exists('sim/true_traj.pkl'):
    TT=solve_ivp(lambda z,x: v_true(np.atleast_1d(x),z)/C.vz_of(z),[zs0,L],x0,t_eval=zeval,rtol=1e-8,atol=1e-11,max_step=0.005).y
    pickle.dump(TT,open('sim/true_traj.pkl','wb'))
TT=pickle.load(open('sim/true_traj.pkl','rb'))

fn='sim/stage2_results.pkl'
res=pickle.load(open(fn,'rb')) if os.path.exists(fn) else {}
if 'clean' not in res:
    Vc=np.array([vfield(mp,None) for mp in MAPS]); Tc=traj_recon(Vc)
    res['clean']=dict(V=Vc.astype(np.float32),rms=np.sqrt(np.mean((Tc-TT)**2,axis=0)),env=Tc[:,-1].std()/TT[:,-1].std())
    pickle.dump(res,open(fn,'wb'))
    print(f"clean-pipeline: final RMS {res['clean']['rms'][-1]*1e6:.1f} um, env ratio {res['clean']['env']:.3f}")
res.setdefault('runs',[])
target=int(sys.argv[1]) if len(sys.argv)>1 else 30
i=len(res['runs'])
while i<target:
    t0=time.time(); rng=np.random.default_rng(1000+i)
    V=np.array([vfield(mp,rng) for mp in MAPS]); T=traj_recon(V)
    res['runs'].append(dict(V=V.astype(np.float32),rms=np.sqrt(np.mean((T-TT)**2,axis=0)),env=T[:,-1].std()/TT[:,-1].std()))
    i+=1
    if i%5==0 or i==target:
        pickle.dump(res,open(fn,'wb')); print(f"  realization {i}/{target} ({time.time()-t0:.1f}s/run) checkpointed")
rms=np.array([r['rms'] for r in res['runs']]); env=np.array([r['env'] for r in res['runs']])
print(f"\n=== STAGE 2 REFERENCE CELL (Na=1e4, {len(res['runs'])} realizations) ===")
print(f"final trajectory RMS : {rms[:,-1].mean()*1e6:6.1f} +- {rms[:,-1].std()*1e6:.1f} um  ({rms[:,-1].mean()/C.FRINGE*100:.1f}% of fringe)")
print(f"clean-pipeline limit : {res['clean']['rms'][-1]*1e6:6.1f} um  ({res['clean']['rms'][-1]/C.FRINGE*100:.1f}%)")
print(f"envelope-width ratio : {env.mean():.3f} +- {env.std():.3f}  (clean {res['clean']['env']:.3f}, true=1)")
Vs=np.array([r['V'] for r in res['runs']]); k=10; selr=np.abs(xr)<400e-6
print(f"per-point sigma_v (z={zj[k]*100:.0f} cm, |x|<400um): {Vs[:,k,selr].std(axis=0).mean()*1e3:.2f} mm/s")
