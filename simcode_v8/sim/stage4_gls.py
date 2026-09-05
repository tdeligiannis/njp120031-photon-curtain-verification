"""Stage 4/D: GLS spline field estimator on the independent per-x_c moment samples.
Model: mom(x) = rho_w(x) + t * d/dx[rho_w v](x),  v = sum_j c_j B_j (cubic B-splines).
GLS with known heteroscedastic sigma_k; mild 2nd-difference ridge. Gate: noiseless truncation bias."""
import numpy as np, pickle, os, sys
sys.path.insert(0,'sim')
from scipy.interpolate import BSpline
import v8_config as C
d,sig0,L=C.d,C.sig0,C.L; W=C.W
NB=2500; xf_edges=np.linspace(-1.25e-3,1.25e-3,NB+1); xf=0.5*(xf_edges[1:]+xf_edges[:-1])
zj=np.linspace(0.05,0.57,20)
MAPS=pickle.load(open('sim/stage3_maps_eta0.pkl','rb'))
MAPS=[MAPS[k] for k in range(20)]

def plane_design(mp):
    """Design matrix, exact moment, and per-point sigma for one plane."""
    xc=mp['xc']; t=mp['t']; lam_f=C.lam_f(mp['z'])
    span=xc[mp['valid']]; lo,hi=span.min(),span.max()
    dk=max(lam_f/4,3*(xc[1]-xc[0]))
    nk=max(int((hi-lo)/dk)+1,5)
    knots=np.concatenate([[lo]*3,np.linspace(lo,hi,nk),[hi]*3])
    nb=nk+2
    B=BSpline.design_matrix(np.clip(xc,lo,hi),knots,3).toarray()   # (151, nb)
    rw=mp['V']*mp['rho_w']+(1-mp['V'])*mp['rho_w_inc']   # measured (mixture) window density
    drw=np.gradient(rw,xc)
    dB=np.array([BSpline(knots,np.eye(nb)[j],3).derivative()(np.clip(xc,lo,hi)) for j in range(nb)]).T
    A=t*(drw[:,None]*B+rw[:,None]*dB)
    prodb=(mp['a_w']*mp['rho_b']).astype(float); lev=xf[None,:]-xc[:,None]
    win=np.abs(prodb)>1e-3*np.abs(prodb).max(axis=1,keepdims=True)
    mom_exact=((lev*prodb).sum(axis=1)+mp['D_corr'])/W
    Mout=(lev*prodb*(~win)).sum(axis=1)
    return dict(xc=xc,t=t,A=A,B=B,knots=knots,nb=nb,win=win,Mout=Mout,mom_exact=mom_exact,
                lo=lo,hi=hi,rho=mp['rho_b'].astype(float),aw=mp['a_w'].astype(float),
                Dc=mp['D_corr'],rw=rw,valid=mp['valid'],
                rho_inc=mp['rho_b_inc'].astype(float),bg=mp['bg_b'].astype(float),aw_inc=mp['a_w_inc'].astype(float),
                V=mp['V'],Psc=mp['Psc'],SIG=mp['SIG'])

def gls_fit(P,y,sig,lam=1e-3):
    Wd=1.0/sig**2
    AtW=P['A'].T*Wd[None,:]
    M2=AtW@P['A']
    R=np.zeros((P['nb'],P['nb']))
    for i in range(1,P['nb']-1):
        v=np.zeros(P['nb']); v[i-1:i+2]=[1,-2,1]; R+=np.outer(v,v)
    sc=np.trace(M2)/max(np.trace(R),1e-30)
    c=np.linalg.solve(M2+lam*sc*R,AtW@y)
    return P['B']@c

def mom_noisy(P,rng,Na):
    rho=P['rho']
    if '_pm' not in P:                                   # v8 mixture, cached per plane
        rhoV=P['V']*rho+(1-P['V'])*P['rho_inc'][None,:]  # m12 visibility mixture
        base=(1-P['Psc'])*rhoV+P['Psc']*P['bg'][None,:]  # per-plane SE background
        sh=lambda s: np.stack([np.interp(xf,xf+s,base[i],left=0,right=0) for i in range(base.shape[0])])
        pm=(1-C.F_IMP)*base+0.5*C.F_IMP*(sh(+C.DSCR_IMP)+sh(-C.DSCR_IMP))   # f-impurity split
        recw=(1-C.F_IMP)*(1-P['Psc'])*(P['V']*P['aw']*rho+(1-P['V'])*P['aw_inc']*P['rho_inc'][None,:])
        P['_pm']=pm; P['_ahm']=np.where(pm>0,recw/np.maximum(pm,1e-300),0.0)
    pm=P['_pm']; ahm=P['_ahm']
    lev=xf[None,:]-P['xc'][:,None]
    pv=np.concatenate([pm,np.maximum(1-pm.sum(1,keepdims=True),0)],axis=1)
    n=np.stack([rng.multinomial(Na,pv[i]/pv[i].sum())[:NB] for i in range(P['xc'].size)])
    ah=ahm+rng.standard_normal(pm.shape)*P['SIG']/np.sqrt(np.maximum(n,1))
    dat=(lev*(n/Na)*ah*P['win']).sum(axis=1)
    dat/= (1-C.F_IMP)*(1-P['Psc'])   # analyst normalization: known factors
    mom=(dat+P['Mout']+P['Dc'])/W
    var=((lev**2)*rho*P['win']).sum(axis=1)*P['SIG']**2/(Na*W**2)
    var/= ((1-C.F_IMP)*(1-P['Psc']))**2   # noise on the normalized moment
    return mom,np.sqrt(var)

if __name__=="__main__":
    if not os.path.exists('sim/stage4_design.pkl'):
        DES=[plane_design(mp) for mp in MAPS]; pickle.dump(DES,open('sim/stage4_design.pkl','wb'))
    else: DES=pickle.load(open('sim/stage4_design.pkl','rb'))
    print("GLS GATE (noiseless): basis-truncation bias of v_GLS vs pipeline v_nl, core region")
    for k in [0,5,9,14,19]:
        P=DES[k]; y=P['mom_exact']-P['rw']
        vg=gls_fit(P,y,np.ones_like(y))            # unweighted on exact data
        vnl_num=np.zeros(P['xc'].size)
        h=P['xc'][1]-P['xc'][0]; integ=P['mom_exact']-P['rw']
        for i in range(1,P['xc'].size):
            vnl_num[i]=vnl_num[i-2]+h/3*(integ[i-2]+4*integ[i-1]+integ[i]) if i%2==0 else vnl_num[i-1]+h/2*(integ[i-1]+integ[i])
        vnl=vnl_num/(P['t']*np.maximum(P['rw'],1e-300))
        core=P['rw']>=0.05*P['rw'].max()
        dev=np.abs(vg-vnl)[core]; scale=np.median(np.abs(vnl[core]))
        print(f"  z={zj[k]*100:5.1f} cm : median|v_GLS - v_nl|/scale = {np.median(dev)/scale*100:6.2f}%  (dof={P['nb']})")
