"""Stage 4/D production: GLS field maps + trajectories at Na sweep; dipole bias through GLS."""
import numpy as np, pickle, sys, time
sys.path.insert(0,'sim')
import stage4_gls as G
gls_fit=G.gls_fit; mom_noisy=G.mom_noisy
from scipy.interpolate import RectBivariateSpline
from sim_core import make_psi
import v8_config as C
d,sig0,L=C.d,C.sig0,C.L; W=C.W
zj=np.linspace(0.05,0.57,20)
DES=pickle.load(open('sim/stage4_design.pkl','rb'))
psi,_,_=make_psi(d,sig0,tmap=C.t_of_z)
TT=pickle.load(open('sim/true_traj.pkl','rb'))
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
        P=DES[k]
        out[k]=np.interp(xr,P['xc'],rows[k],left=0,right=0)
        out[k][(xr<P['lo'])|(xr>P['hi'])]=0
    return out
# noiseless GLS truth field for reference and error normalization
Vt=[gls_fit(DES[k],DES[k]['mom_exact']-DES[k]['rw'],np.ones(DES[k]['xc'].size)) for k in range(20)]
res={'zj':zj}
# matched-filter detection SNR per plane at Na=1e4 (scales as sqrt(Na))
print("=== matched-filter field-detection SNR per plane (Na=1e4; scale by sqrt(Na/1e4)) ===")
snr=np.zeros(20)
for k in range(20):
    P=DES[k]; _,sg=mom_noisy(P,np.random.default_rng(0),int(1e4))
    snr[k]=np.sqrt((((P['mom_exact']-P['rw'])/sg)**2).sum())
res['snr_1e4']=snr
print("  ".join(f"z{zj[k]*100:.0f}:{snr[k]:.2f}" for k in range(0,20,3)))
print(f"  global (all planes pooled): {np.sqrt((snr**2).sum()):.2f}")
# Na sweep: 5 realizations each
print("\n=== GLS field + trajectories vs Na (5 realizations) ===")
res['sweep']={}
for Na in [1e4,1e5,1e6,1e7]:
    t0=time.time(); fr=[]; tr=[]
    for r in range(5):
        rng=np.random.default_rng(4000+r)
        rows=[]
        ferr=[]
        for k in range(20):
            P=DES[k]; mom,sg=mom_noisy(P,rng,int(Na))
            vg=gls_fit(P,mom-P['rw'],sg)
            rows.append(vg)
            core=P['rw']>=0.05*P['rw'].max()
            ferr.append(np.sqrt(np.mean((vg-Vt[k])[core]**2))/np.median(np.abs(Vt[k][core])))
        fr.append(ferr)
        T=traj(resample(rows)); tr.append(np.sqrt(np.mean((T-TT)**2,axis=0))[-1])
    fr=np.array(fr); tr=np.array(tr)
    res['sweep'][Na]=dict(field_rel=fr,traj=tr)
    print(f"Na={Na:.0e}: traj RMS {np.mean(tr)*1e6:8.1f}+-{np.std(tr)*1e6:6.1f} um ({np.mean(tr)/C.FRINGE*100:6.1f}% fringe) | field rel-err: z=5cm {fr[:,0].mean()*100:5.1f}%  z=30cm {fr[:,9].mean()*100:6.1f}%  z=57cm {fr[:,19].mean()*100:6.1f}%  [{time.time()-t0:.0f}s]")
# dipole bias through GLS (linear in mom): momD per cell -> GLS -> trajectories
print("\n=== dipole residual bias through GLS pipeline ===")
def momex(nm):
    Mm=pickle.load(open(f'sim/stage3_maps_{nm}.pkl','rb')); Mm=[Mm[k] for k in range(20)]
    NBl=2500; xfl=G.xf
    out=[]
    for k in range(20):
        mp=Mm[k]; prodb=(mp['a_w']*mp['rho_b']).astype(float)
        lev=xfl[None,:]-mp['xc'][:,None]
        out.append(((lev*prodb).sum(axis=1)+mp['D_corr'])/W)
    return out
mE={nm:momex(nm) for nm in ['eta1em3','eta1em2','eta1em3off','eta1em2off']}
Bfield={}
V0=resample(Vt)
T0=traj(V0)
res['bias_gls']={}
for eb,er,Dm in [('1e-2','1e-3',[mE['eta1em2'][k]-mE['eta1em2off'][k] for k in range(20)]),
                 ('1e-3','1e-3',[mE['eta1em3'][k]-mE['eta1em3off'][k] for k in range(20)]),
                 ('1e-2','1e-4',[0.1*(mE['eta1em2'][k]-mE['eta1em2off'][k]) for k in range(20)]),
                 ('1e-3','1e-4',[0.1*(mE['eta1em3'][k]-mE['eta1em3off'][k]) for k in range(20)])]:
    rows=[Vt[k]+gls_fit(DES[k],Dm[k],np.ones(DES[k]['xc'].size)) for k in range(20)]
    T=traj(resample(rows)); r=np.sqrt(np.mean((T-T0)**2,axis=0))[-1]
    res['bias_gls'][(eb,er)]=r
    print(f"  eta={eb}, eps={er} (prod {float(eb)*float(er):.0e}): GLS-pipeline traj bias {r*1e6:8.2f} um ({r/C.FRINGE*100:.3f}% fringe)")
pickle.dump(res,open('sim/stage4_results.pkl','wb')); print("\nsaved sim/stage4_results.pkl")
