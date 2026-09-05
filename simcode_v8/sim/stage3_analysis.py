"""Stage 3 step B: deterministic dipole-residual bias fields, budget validation, product test."""
import numpy as np, pickle, sys
sys.path.insert(0,'sim')
import v8_config as C
W=C.W
NB=2500; xf_edges=np.linspace(-1.25e-3,1.25e-3,NB+1); xf=0.5*(xf_edges[1:]+xf_edges[:-1])
zj=np.linspace(0.05,0.57,20)
M={nm:pickle.load(open(f'sim/stage3_maps_{nm}.pkl','rb')) for nm in ['eta0','eta1em3','eta1em2','eta1em3off','eta1em2off']}
from sim_core import cum_integrate
def v_nl(mp):
    prodb=mp['a_w']*mp['rho_b']; lev=xf[None,:]-mp['xc'][:,None]
    win=np.abs(prodb)>1e-3*np.abs(prodb).max(axis=1,keepdims=True)
    mom=((lev*prodb*win).sum(axis=1)+(lev*prodb*(~win)).sum(axis=1)+mp['D_corr'])/W
    v=cum_integrate(mom-mp['rho_w'],mp['xc'])/(mp['t']*np.maximum(mp['rho_w'],1e-300))
    v[~mp['valid']]=0.0
    return v
V={nm:np.array([v_nl(M[nm][k]) for k in range(20)]) for nm in M}
# per-plane sigma_v from stage2 (noise floor at Na=1e4), rho-weighted over valid support
S2=pickle.load(open('sim/stage2_results.pkl','rb'))
xr=np.linspace(-950e-6,950e-6,381)
Vs=np.array([r['V'] for r in S2['runs']])
sigv=np.zeros(20)
for k in range(20):
    mp=M['eta0'][k]; sel=np.interp(xr,mp['xc'],mp['rho_w'])>1e-2*mp['rho_w'].max()
    sigv[k]=Vs[:,k,sel].std(axis=0).mean()
res={'zj':zj,'sigv_1e4':sigv}
print("=== raw wiggle W_eta = v_nl[eta]-v_nl[0]  and residual D_eta = v_nl[eta]-v_nl[0.999 eta] ===")
print("(rho_w-weighted RMS over valid, um/s)")
print(f"{'z(cm)':>6} {'wig 1e-3':>9} {'wig 1e-2':>9} {'D 1e-3':>8} {'D 1e-2':>8} {'wigscale':>9} {'sig_v(1e4)':>11}")
tab=np.zeros((20,4))
for k in range(20):
    mp=M['eta0'][k]; wt=mp['rho_w']*mp['valid']; wt/=wt.sum()
    rms=lambda A: np.sqrt((wt*A[k]**2).sum())
    W3=V['eta1em3']-V['eta0']; W2=V['eta1em2']-V['eta0']
    D3=V['eta1em3']-V['eta1em3off']; D2=V['eta1em2']-V['eta1em2off']
    tab[k]=[rms(W3),rms(W2),rms(D3),rms(D2)]
    if k in [0,4,9,14,19]:
        print(f"{zj[k]*100:6.1f} {tab[k,0]*1e6:9.2f} {tab[k,1]*1e6:9.2f} {tab[k,2]*1e6:8.3f} {tab[k,3]*1e6:8.3f} {C.d/C.t_of_z(zj[k])*1e6:9.2f} {sigv[k]*1e6:11.0f}")
res['wiggle_1em3'],res['wiggle_1em2'],res['D_1em3'],res['D_1em2']=tab.T
# product-only test: (eta=1e-3,eps=1e-3) vs (eta=1e-2,eps=1e-4), both product 1e-6
print("\n=== product-only test: D_{1e-3} vs 0.1*D_{1e-2} (both eta*eps=1e-6) ===")
for k in [0,9,19]:
    mp=M['eta0'][k]; wt=mp['rho_w']*mp['valid']; wt/=wt.sum()
    D3=(V['eta1em3']-V['eta1em3off'])[k]; D2s=0.1*(V['eta1em2']-V['eta1em2off'])[k]
    num=np.sqrt((wt*(D3-D2s)**2).sum()); den=np.sqrt((wt*D3**2).sum())
    print(f"  z={zj[k]*100:5.1f} cm : RMS(diff)/RMS(D_1e-3) = {num/den*100:6.2f}%")
# cell table: residuals vs budget (10% of wiggle) and vs noise
print("\n=== cell residual summary (worst plane, rho-weighted RMS) ===")
cells=[('1e-2','1e-3',res['D_1em2']),( '1e-2','1e-4',0.1*res['D_1em2']),('1e-3','1e-3',res['D_1em3']),('1e-3','1e-4',0.1*res['D_1em3'])]
wig=np.array([C.d/C.t_of_z(z) for z in zj])   # v8 wiggle scale d/t(z) per plane
for eb,er,Dv in cells:
    frac=(Dv/wig); kw=np.argmax(frac)
    print(f"  eta={eb}, eps={er} (prod {float(eb)*float(er):.0e}): worst D/wigscale = {frac[kw]*100:6.2f}% at z={zj[kw]*100:.0f}cm ; D/sigma_v(1e4) worst = {(Dv/sigv).max()*100:.3f}%")
pickle.dump(res,open('sim/stage3_bias_results.pkl','wb')); print("\nsaved sim/stage3_bias_results.pkl")
