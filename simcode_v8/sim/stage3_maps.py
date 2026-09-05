"""Stage 3 maps: dipole-residual phase masks on a common fine grid.
Masks: eta in {0, 1e-3, 1e-2} (truth) + {0.999e-3, 0.999e-2} (epsilon-derivative offsets).
Grid: +-6mm, 2^17 (dx=0.0916um; Nyquist 3.4e7 > mask gradient 1.9e7 at eta=1e-2;
kicked-stripe max displacement 5.6mm inside grid, wrap tail <1e-4 of stripe mass).
Usage: python3 sim/stage3_maps.py <mask_index> <plane_lo> <plane_hi>"""
import numpy as np, pickle, os, sys, time
sys.path.insert(0,'sim'); from sim_core import make_psi, make_psi_parts, hbar, M
import v8_config as C
from scipy.signal import fftconvolve

d,sig0,L = C.d,C.sig0,C.L; W=C.W; MASK_REL=1e-4   # v8: per-plane GA=ga_m(z) inside maps_for_plane
ETAS=[0.0,1e-3,1e-2,0.999e-3,0.999e-2]
MNAMES=['eta0','eta1em3','eta1em2','eta1em3off','eta1em2off']
zj=np.linspace(0.05,0.57,20)
NB=2500; xf_edges=np.linspace(-1.25e-3,1.25e-3,NB+1); xf=0.5*(xf_edges[1:]+xf_edges[:-1])
Xh,N=6e-3,2**17; xg=np.linspace(-Xh,Xh,N,endpoint=False); dx=xg[1]-xg[0]
psi,dpsi,v_true=make_psi(d,sig0,tmap=C.t_of_z)
psa,psb=make_psi_parts(d,sig0,tmap=C.t_of_z)
kk=2*np.pi*np.fft.fftfreq(N,dx)
ilo=np.searchsorted(xg,xf_edges[:-1]); ihi=np.searchsorted(xg,xf_edges[1:])
REF=pickle.load(open('sim/stage2_maps.pkl','rb'))     # reuse per-plane adaptive xc grids
def binsum2(arr):
    cs=np.concatenate([np.zeros(arr.shape[:-1]+(1,)),np.cumsum(arr,axis=-1)],axis=-1)
    return cs[...,ihi]-cs[...,ilo]

def maps_for_plane(z, eta, xc):
    t=C.T-C.t_of_z(z); GA=C.ga_m_of(z); kern=np.exp(-1j*hbar*kk**2*t/(2*M))
    psr=psi(xg,z); nrm=np.sqrt(np.sum(np.abs(psr)**2)*dx); ps=psr/nrm
    pa=psa(xg,z)/nrm; pb=psb(xg,z)/nrm
    rho_j=np.abs(ps)**2
    psL0=np.fft.ifft(np.fft.fft(ps)*kern); rho_L0=np.abs(psL0)**2
    paL=np.fft.ifft(np.fft.fft(pa)*kern); pbL=np.fft.ifft(np.fft.fft(pb)*kern)
    rho_inc=np.abs(paL)**2+np.abs(pbL)**2; rho_inc/=np.sum(rho_inc)*dx
    Dk=C.vrec*(C.T-C.t_of_z(z)); nb_box=max(int(2*Dk/dx)|1,3)
    bg=fftconvolve(rho_L0,np.ones(nb_box)/nb_box,mode='same')
    rw=np.array([np.sum(np.exp(-(xg-c)**2/W**2)*rho_j)*dx for c in xc])   # mask leaves |psi(z_j)| unchanged
    ctrg=xf[np.clip(np.searchsorted(xf_edges,xg,side='right')-1,0,NB-1)]  # side='right': match binsum's left-inclusive bins (x=0 sits exactly on an edge)
    a_w=np.zeros((xc.size,NB),np.float32); rho_b=np.zeros((xc.size,NB),np.float32)
    D_corr=np.zeros(xc.size); numd=np.zeros(xc.size); p_out=np.zeros(xc.size)
    CH=40
    for s0 in range(0,xc.size,CH):
        cs=xc[s0:s0+CH]
        f=np.exp(-(xg[None,:]-cs[:,None])**2/W**2)
        mask=np.exp(1j*eta*GA*f) if eta!=0 else 1.0
        base=(f*ps[None,:]) if eta==0 else (f*(mask*ps[None,:]))
        psLp=np.fft.ifft(np.fft.fft((mask*ps[None,:]) if eta!=0 else np.broadcast_to(ps,(cs.size,N)),axis=1)*kern[None,:],axis=1)
        G =np.fft.ifft(np.fft.fft(base,axis=1)*kern[None,:],axis=1)
        Gd=np.fft.ifft(np.fft.fft((2*(xg[None,:]-cs[:,None])/W**2)*base,axis=1)*kern[None,:],axis=1)
        rhoL=np.abs(psLp)**2
        rho_b[s0:s0+CH]=binsum2(rhoL*dx)
        prodfull=W*np.real(Gd*np.conj(psLp))*dx
        pb=binsum2(prodfull)
        a_w[s0:s0+CH]=np.where(rho_b[s0:s0+CH]>0,pb/np.maximum(rho_b[s0:s0+CH],1e-30),0.0)  # 1e-30: representable in float32
        D_corr[s0:s0+CH]=((xg-ctrg)[None,:]*prodfull).sum(axis=1)
        numd[s0:s0+CH]=((xg[None,:]-cs[:,None])*np.real(np.conj(psLp)*G)).sum(axis=1)*dx
        p_out[s0:s0+CH]=1-rho_b[s0:s0+CH].sum(axis=1)
    vfor=np.where(rw>0,numd/(t*np.maximum(rw,1e-300)),np.nan)
    return dict(z=z,t=t,xc=xc,a_w=a_w,rho_b=rho_b,rho_w=rw,D_corr=D_corr,num_direct=numd,
                v_formula=vfor,valid=rw>=MASK_REL*rw.max(),p_out=p_out,
                rho_b_inc=binsum2(rho_inc*dx),bg_b=binsum2(bg*dx),
                V=float(C.V_of(z)),Psc=float(C.Psc_of(z)),SIG=float(C.SIG_of(z)))

if __name__=="__main__":
    mi,lo,hi=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])
    eta=ETAS[mi]; nm=MNAMES[mi]
    fn=f'sim/stage3_maps_{nm}.pkl'
    store=pickle.load(open(fn,'rb')) if os.path.exists(fn) else {}
    for k in range(lo,hi):
        if k in store: continue
        t0=time.time(); store[k]=maps_for_plane(zj[k],eta,REF[k]['xc'])
        store[k]['a_w_inc']=REF[k]['a_w_inc']; store[k]['rho_w_inc']=REF[k]['rho_w_inc']  # eta-independent at (1-V)*eta order; same xc grid
        pickle.dump(store,open(fn,'wb'))
        print(f"[{nm}] plane {k} (z={zj[k]*100:.1f}cm) {time.time()-t0:.1f}s  max p_out={store[k]['p_out'].max()*100:.3f}%")
    if len(store)==20: print(f"[{nm}] COMPLETE")
