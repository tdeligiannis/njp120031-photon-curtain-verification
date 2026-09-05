"""Stage 2 maps v3 — adaptive per-plane x_c grids (151 pts spanning the rho_w support), 1um xf bins."""
import numpy as np, pickle, os, sys, time
sys.path.insert(0,'sim'); from sim_core import make_psi, make_psi_parts, hbar, M
import v8_config as C
from scipy.signal import fftconvolve

d,sig0,L = C.d,C.sig0,C.L; W=C.W   # v8: transverse kernel width = w_eff (convention D3)
MASK_REL=1e-4; GRID_REL=1e-8; NXC=151   # GRID_REL: integration-domain tail (cumulative needs N(xc[0])~0); MASK_REL: analysis validity
zj=np.linspace(0.05,0.57,20)
NB=2500; xf_edges=np.linspace(-1.25e-3,1.25e-3,NB+1); xf=0.5*(xf_edges[1:]+xf_edges[:-1])
Xh,N=3e-3,2**15; xg=np.linspace(-Xh,Xh,N,endpoint=False); dx=xg[1]-xg[0]
psi,dpsi,v_true=make_psi(d,sig0,tmap=C.t_of_z)
psa,psb=make_psi_parts(d,sig0,tmap=C.t_of_z)
kk=2*np.pi*np.fft.fftfreq(N,dx)
ilo=np.searchsorted(xg,xf_edges[:-1]); ihi=np.searchsorted(xg,xf_edges[1:])
def binsum(arr):
    cs=np.concatenate([np.zeros(arr.shape[:-1]+(1,)),np.cumsum(arr,axis=-1)],axis=-1)
    return cs[...,ihi]-cs[...,ilo]

def maps_for_plane(z):
    t=C.T-C.t_of_z(z); kern=np.exp(-1j*hbar*kk**2*t/(2*M))
    psr=psi(xg,z); nrm=np.sqrt(np.sum(np.abs(psr)**2)*dx); ps=psr/nrm
    pa=psa(xg,z)/nrm; pb=psb(xg,z)/nrm; pa_s,pb_s=pa,pb   # capture: pb is reused below for the binned product
    rho_j=np.abs(ps)**2
    # adaptive x_c grid: span of w-blurred support above MASK_REL, padded 3w
    rw_fine=np.convolve(rho_j,np.exp(-(np.arange(-8e-6,8e-6,dx))**2/W**2),mode='same')*dx
    sup=xg[rw_fine>=GRID_REL*rw_fine.max()]
    xc=np.linspace(sup.min()-3*C.w_spot,sup.max()+3*C.w_spot,NXC)
    psL=np.fft.ifft(np.fft.fft(ps)*kern); rho_L=np.abs(psL)**2
    paL=np.fft.ifft(np.fft.fft(pa)*kern); pbL=np.fft.ifft(np.fft.fft(pb)*kern)
    rho_inc=np.abs(paL)**2+np.abs(pbL)**2; rho_inc/=np.sum(rho_inc)*dx      # m12 incoherent part
    Dk=C.vrec*(C.T-C.t_of_z(z)); nb_box=max(int(2*Dk/dx)|1,3)
    bg=fftconvolve(rho_L,np.ones(nb_box)/nb_box,mode='same')                # SE-recoil-smeared bg
    rho_b=binsum(rho_L*dx); rho_b_inc=binsum(rho_inc*dx); bg_b=binsum(bg*dx)
    Fm=np.exp(-(xg[None,:]-xc[:,None])**2/W**2)
    G =np.fft.ifft(np.fft.fft(Fm*ps[None,:],axis=1)*kern[None,:],axis=1)
    Gd=np.fft.ifft(np.fft.fft((2*(xg[None,:]-xc[:,None])/W**2)*Fm*ps[None,:],axis=1)*kern[None,:],axis=1)
    rw=(Fm*rho_j[None,:]).sum(axis=1)*dx
    rho_j_inc=np.abs(pa_s)**2+np.abs(pb_s)**2; rho_j_inc/=np.sum(rho_j_inc)*dx
    rw_inc=(Fm*rho_j_inc[None,:]).sum(axis=1)*dx   # incoherent window density at the plane
    numd=((xg[None,:]-xc[:,None])*np.real(np.conj(psL)[None,:]*G)).sum(axis=1)*dx
    prodfull=W*np.real(Gd*np.conj(psL)[None,:])*dx
    pb=binsum(prodfull)
    # m12 PROPER incoherent record model: per-slit products (no cross terms).
    # Using the coherent a_w on rho_inc is wrong at fringe minima (weak-value peaks x envelope).
    Gda=np.fft.ifft(np.fft.fft((2*(xg[None,:]-xc[:,None])/W**2)*Fm*pa_s[None,:],axis=1)*kern[None,:],axis=1)
    prod_inc=W*np.real(Gda*np.conj(paL)[None,:])*dx
    Gda=np.fft.ifft(np.fft.fft((2*(xg[None,:]-xc[:,None])/W**2)*Fm*pb_s[None,:],axis=1)*kern[None,:],axis=1)
    prod_inc=prod_inc+W*np.real(Gda*np.conj(pbL)[None,:])*dx; del Gda
    prod_inc/= np.sum(np.abs(paL)**2+np.abs(pbL)**2)*dx   # same rescale as rho_inc normalization
    pb_inc=binsum(prod_inc)
    ctrg=xf[np.clip(np.searchsorted(xf_edges,xg,side='right')-1,0,NB-1)]  # side='right': match binsum's left-inclusive bins (x=0 sits exactly on an edge)
    D_corr=((xg-ctrg)[None,:]*prodfull).sum(axis=1)
    a_w=np.where(rho_b[None,:]>0,pb/np.maximum(rho_b[None,:],1e-300),0.0)
    a_w_inc=np.where(rho_b_inc[None,:]>0,pb_inc/np.maximum(rho_b_inc[None,:],1e-300),0.0)
    vfor=np.where(rw>0,numd/(t*np.maximum(rw,1e-300)),np.nan)
    return dict(z=z,t=t,xc=xc,a_w=a_w.astype(np.float32),rho_b=rho_b,rho_w=rw,v_formula=vfor,
                num_direct=numd,D_corr=D_corr,v_true=v_true(xc,z),valid=rw>=MASK_REL*rw.max(),p_out=1-rho_b.sum(),
                rho_b_inc=rho_b_inc,bg_b=bg_b,a_w_inc=a_w_inc.astype(np.float32),rho_w_inc=rw_inc,
                V=float(C.V_of(z)),Psc=float(C.Psc_of(z)),SIG=float(C.SIG_of(z)))

if __name__=="__main__":
    lo,hi=int(sys.argv[1]),int(sys.argv[2])
    fn='sim/stage2_maps_partial.pkl'
    store=pickle.load(open(fn,'rb')) if os.path.exists(fn) else {}
    for k in range(lo,hi):
        if k in store: continue
        t0=time.time(); store[k]=maps_for_plane(zj[k]); pickle.dump(store,open(fn,'wb'))
        print(f"plane {k} (z={zj[k]*100:.1f} cm) dxc={(store[k]['xc'][1]-store[k]['xc'][0])*1e6:.2f} um done {time.time()-t0:.1f}s")
    if len(store)==20:
        pickle.dump([store[k] for k in range(20)],open('sim/stage2_maps.pkl','wb')); print("FULL MAP SET SAVED")
