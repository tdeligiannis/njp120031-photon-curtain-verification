"""Build simcode_v8/ from the v7 bundle: install v8_config.py + rewritten sim_core.py +
gate-only stage2_reference.py, patch the eight stage scripts (asserted-unique edits),
purge stale v7 pickles. Every edit tagged; any count!=1 aborts the build."""
import os, shutil, glob

SRC = "/home/claude/njp_v8/simcode"
DST = "/home/claude/njp_v8/simcode_v8"
NEW = "/home/claude/njp_v8/simcode_v8_new"

if os.path.exists(DST): shutil.rmtree(DST)
shutil.copytree(SRC, DST)
os.makedirs(f"{DST}/legacy", exist_ok=True)
for p in glob.glob(f"{DST}/sim/*.pkl"): os.remove(p)
shutil.move(f"{DST}/NJP120031_final_verification_output.txt",
            f"{DST}/legacy/NJP120031_v7_final_verification_output.txt")
shutil.copy(f"{NEW}/v8_config.py", f"{DST}/sim/v8_config.py")
shutil.copy(f"{NEW}/sim_core.py",  f"{DST}/sim/sim_core.py")
print("copied tree, purged v7 pickles, installed v8_config + sim_core")

def patch(fn, edits):
    p = f"{DST}/sim/{fn}"
    s = open(p).read()
    for tag, old, new in edits:
        n = s.count(old)
        assert n == 1, f"[{fn}:{tag}] count={n}"
        s = s.replace(old, new)
        print(f"  {fn}:{tag} OK")
    open(p, "w").write(s)

# ============================ stage2_maps.py ============================
patch("stage2_maps.py", [
("imports",
 "sys.path.insert(0,'sim'); from sim_core import make_psi, hbar, M",
 "sys.path.insert(0,'sim'); from sim_core import make_psi, make_psi_parts, hbar, M\nimport v8_config as C\nfrom scipy.signal import fftconvolve"),
("params",
 "d,sig0,w,vz,L = 2.5e-6,0.6e-6,2e-6,1.5,0.60",
 "d,sig0,L = C.d,C.sig0,C.L; W=C.W   # v8: transverse kernel width = w_eff (convention D3)"),
("psi",
 "psi,dpsi,v_true=make_psi(d,sig0,vz)",
 "psi,dpsi,v_true=make_psi(d,sig0,tmap=C.t_of_z)\npsa,psb=make_psi_parts(d,sig0,tmap=C.t_of_z)"),
("tmap",
 "def maps_for_plane(z):\n    t=(L-z)/vz; kern=np.exp(-1j*hbar*kk**2*t/(2*M))",
 "def maps_for_plane(z):\n    t=C.T-C.t_of_z(z); kern=np.exp(-1j*hbar*kk**2*t/(2*M))"),
("perslit-norm",
 "    ps=psi(xg,z); ps=ps/np.sqrt(np.sum(np.abs(ps)**2)*dx)\n    rho_j=np.abs(ps)**2",
 "    psr=psi(xg,z); nrm=np.sqrt(np.sum(np.abs(psr)**2)*dx); ps=psr/nrm\n"
 "    pa=psa(xg,z)/nrm; pb=psb(xg,z)/nrm\n    rho_j=np.abs(ps)**2"),
("adaptive-kernel",
 "rw_fine=np.convolve(rho_j,np.exp(-(np.arange(-8e-6,8e-6,dx))**2/w**2),mode='same')*dx",
 "rw_fine=np.convolve(rho_j,np.exp(-(np.arange(-8e-6,8e-6,dx))**2/W**2),mode='same')*dx"),
("pad",
 "xc=np.linspace(sup.min()-3*w,sup.max()+3*w,NXC)",
 "xc=np.linspace(sup.min()-3*C.w_spot,sup.max()+3*C.w_spot,NXC)"),
("screen-extras",
 "    psL=np.fft.ifft(np.fft.fft(ps)*kern); rho_L=np.abs(psL)**2\n    rho_b=binsum(rho_L*dx)",
 "    psL=np.fft.ifft(np.fft.fft(ps)*kern); rho_L=np.abs(psL)**2\n"
 "    paL=np.fft.ifft(np.fft.fft(pa)*kern); pbL=np.fft.ifft(np.fft.fft(pb)*kern)\n"
 "    rho_inc=np.abs(paL)**2+np.abs(pbL)**2; rho_inc/=np.sum(rho_inc)*dx      # m12 incoherent part\n"
 "    Dk=C.vrec*(C.T-C.t_of_z(z)); nb_box=max(int(2*Dk/dx)|1,3)\n"
 "    bg=fftconvolve(rho_L,np.ones(nb_box)/nb_box,mode='same')                # SE-recoil-smeared bg\n"
 "    rho_b=binsum(rho_L*dx); rho_b_inc=binsum(rho_inc*dx); bg_b=binsum(bg*dx)"),
("window",
 "Fm=np.exp(-(xg[None,:]-xc[:,None])**2/w**2)",
 "Fm=np.exp(-(xg[None,:]-xc[:,None])**2/W**2)"),
("deriv",
 "Gd=np.fft.ifft(np.fft.fft((2*(xg[None,:]-xc[:,None])/w**2)*Fm*ps[None,:],axis=1)*kern[None,:],axis=1)",
 "Gd=np.fft.ifft(np.fft.fft((2*(xg[None,:]-xc[:,None])/W**2)*Fm*ps[None,:],axis=1)*kern[None,:],axis=1)"),
("prod",
 "prodfull=w*np.real(Gd*np.conj(psL)[None,:])*dx",
 "prodfull=W*np.real(Gd*np.conj(psL)[None,:])*dx"),
("dict",
 "    return dict(z=z,t=t,xc=xc,a_w=a_w.astype(np.float32),rho_b=rho_b,rho_w=rw,v_formula=vfor,\n"
 "                num_direct=numd,D_corr=D_corr,v_true=v_true(xc,z),valid=rw>=MASK_REL*rw.max(),p_out=1-rho_b.sum())",
 "    return dict(z=z,t=t,xc=xc,a_w=a_w.astype(np.float32),rho_b=rho_b,rho_w=rw,v_formula=vfor,\n"
 "                num_direct=numd,D_corr=D_corr,v_true=v_true(xc,z),valid=rw>=MASK_REL*rw.max(),p_out=1-rho_b.sum(),\n"
 "                rho_b_inc=rho_b_inc,bg_b=bg_b,V=float(C.V_of(z)),Psc=float(C.Psc_of(z)),SIG=float(C.SIG_of(z)))"),
])

# ============================ stage2_run.py ============================
patch("stage2_run.py", [
("imports",
 "sys.path.insert(0,'sim'); from sim_core import make_psi, hbar, M",
 "sys.path.insert(0,'sim'); from sim_core import make_psi, hbar, M\nimport v8_config as C"),
("params",
 "d,sig0,w,vz,L = 2.5e-6,0.6e-6,2e-6,1.5,0.60\n"
 "ga_m,Ng_m,ga_c,Ng_c = 4402.0,1.684e9,-4402.0,4.499e9\n"
 "sig_arm=lambda Ng,ga: np.sqrt(Ng)/(np.sqrt(2)*abs(ga))\n"
 "SIG=1.0/np.sqrt(sig_arm(Ng_m,ga_m)**-2+sig_arm(Ng_c,ga_c)**-2)\n"
 "F_SE=0.005; Na=int(1e4)",
 "d,sig0,L = C.d,C.sig0,C.L; W=C.W\n"
 "Na=int(1e4)   # v8: per-plane SIG/V/Psc carried in the maps; SE bg = Psc(z) x recoil-smeared rho"),
("psi",
 "psi,dpsi,v_true=make_psi(d,sig0,vz)",
 "psi,dpsi,v_true=make_psi(d,sig0,tmap=C.t_of_z)"),
("uw",
 "UW=np.where((xf>=-0.5e-3)&(xf<=0.5e-3),1.0,0.0); UW/=UW.sum()\n",
 ""),
("mixture",
 "    rho=mp['rho_b']; pmix=(1-F_SE)*rho+F_SE*rho.sum()*UW\n"
 "    sigfrac=np.where(pmix>0,(1-F_SE)*rho/np.maximum(pmix,1e-300),0.0)",
 "    rhoV=mp['V']*mp['rho_b']+(1-mp['V'])*mp['rho_b_inc']            # m12 visibility mixture\n"
 "    base=(1-mp['Psc'])*rhoV+mp['Psc']*mp['bg_b']                    # per-plane SE background\n"
 "    sh=lambda s: np.interp(xf,xf+s,base,left=0,right=0)\n"
 "    pmix=(1-C.F_IMP)*base+0.5*C.F_IMP*(sh(+C.DSCR_IMP)+sh(-C.DSCR_IMP))   # f-impurity split\n"
 "    sigfrac=np.where(pmix>0,(1-C.F_IMP)*(1-mp['Psc'])*rhoV/np.maximum(pmix,1e-300),0.0)"),
("sig",
 "ah=sigfrac[None,:]*mp['a_w']+rng.standard_normal((mp['xc'].size,NB))*SIG/np.sqrt(np.maximum(n,1))",
 "ah=sigfrac[None,:]*mp['a_w']+rng.standard_normal((mp['xc'].size,NB))*mp['SIG']/np.sqrt(np.maximum(n,1))"),
("mom",
 "mom=(dat+Mout+mp['D_corr'])/w",
 "mom=(dat+Mout+mp['D_corr'])/W"),
("traj",
 "f=lambda z,x: sp(z,np.clip(x,xr[0],xr[-1]),grid=False)/vz",
 "f=lambda z,x: sp(z,np.clip(x,xr[0],xr[-1]),grid=False)/C.vz_of(z)"),
("truetraj",
 "TT=solve_ivp(lambda z,x: v_true(np.atleast_1d(x),z)/vz,[zs0,L],x0,t_eval=zeval,rtol=1e-8,atol=1e-11,max_step=0.005).y",
 "TT=solve_ivp(lambda z,x: v_true(np.atleast_1d(x),z)/C.vz_of(z),[zs0,L],x0,t_eval=zeval,rtol=1e-8,atol=1e-11,max_step=0.005).y"),
("fringe1",
 "print(f\"final trajectory RMS : {rms[:,-1].mean()*1e6:6.1f} +- {rms[:,-1].std()*1e6:.1f} um  ({rms[:,-1].mean()/735e-6*100:.1f}% of fringe)\")",
 "print(f\"final trajectory RMS : {rms[:,-1].mean()*1e6:6.1f} +- {rms[:,-1].std()*1e6:.1f} um  ({rms[:,-1].mean()/C.FRINGE*100:.1f}% of fringe)\")"),
("fringe2",
 "print(f\"clean-pipeline limit : {res['clean']['rms'][-1]*1e6:6.1f} um  ({res['clean']['rms'][-1]/735e-6*100:.1f}%)\")",
 "print(f\"clean-pipeline limit : {res['clean']['rms'][-1]*1e6:6.1f} um  ({res['clean']['rms'][-1]/C.FRINGE*100:.1f}%)\")"),
])

# ============================ stage3_maps.py ============================
patch("stage3_maps.py", [
("imports",
 "sys.path.insert(0,'sim'); from sim_core import make_psi, hbar, M",
 "sys.path.insert(0,'sim'); from sim_core import make_psi, make_psi_parts, hbar, M\nimport v8_config as C\nfrom scipy.signal import fftconvolve"),
("params",
 "d,sig0,w,vz,L = 2.5e-6,0.6e-6,2e-6,1.5,0.60\nGA=4402.0; MASK_REL=1e-4",
 "d,sig0,L = C.d,C.sig0,C.L; W=C.W; MASK_REL=1e-4   # v8: per-plane GA=ga_m(z) inside maps_for_plane"),
("psi",
 "psi,dpsi,v_true=make_psi(d,sig0,vz)",
 "psi,dpsi,v_true=make_psi(d,sig0,tmap=C.t_of_z)\npsa,psb=make_psi_parts(d,sig0,tmap=C.t_of_z)"),
("tmap",
 "def maps_for_plane(z, eta, xc):\n    t=(L-z)/vz; kern=np.exp(-1j*hbar*kk**2*t/(2*M))",
 "def maps_for_plane(z, eta, xc):\n    t=C.T-C.t_of_z(z); GA=C.ga_m_of(z); kern=np.exp(-1j*hbar*kk**2*t/(2*M))"),
("perslit",
 "    ps=psi(xg,z); ps=ps/np.sqrt(np.sum(np.abs(ps)**2)*dx)\n    rho_j=np.abs(ps)**2",
 "    psr=psi(xg,z); nrm=np.sqrt(np.sum(np.abs(psr)**2)*dx); ps=psr/nrm\n"
 "    pa=psa(xg,z)/nrm; pb=psb(xg,z)/nrm\n    rho_j=np.abs(ps)**2\n"
 "    psL0=np.fft.ifft(np.fft.fft(ps)*kern); rho_L0=np.abs(psL0)**2\n"
 "    paL=np.fft.ifft(np.fft.fft(pa)*kern); pbL=np.fft.ifft(np.fft.fft(pb)*kern)\n"
 "    rho_inc=np.abs(paL)**2+np.abs(pbL)**2; rho_inc/=np.sum(rho_inc)*dx\n"
 "    Dk=C.vrec*(C.T-C.t_of_z(z)); nb_box=max(int(2*Dk/dx)|1,3)\n"
 "    bg=fftconvolve(rho_L0,np.ones(nb_box)/nb_box,mode='same')"),
("rw",
 "rw=np.array([np.sum(np.exp(-(xg-c)**2/w**2)*rho_j)*dx for c in xc])",
 "rw=np.array([np.sum(np.exp(-(xg-c)**2/W**2)*rho_j)*dx for c in xc])"),
("window",
 "        f=np.exp(-(xg[None,:]-cs[:,None])**2/w**2)",
 "        f=np.exp(-(xg[None,:]-cs[:,None])**2/W**2)"),
("deriv",
 "Gd=np.fft.ifft(np.fft.fft((2*(xg[None,:]-cs[:,None])/w**2)*base,axis=1)*kern[None,:],axis=1)",
 "Gd=np.fft.ifft(np.fft.fft((2*(xg[None,:]-cs[:,None])/W**2)*base,axis=1)*kern[None,:],axis=1)"),
("prod",
 "prodfull=w*np.real(Gd*np.conj(psLp))*dx",
 "prodfull=W*np.real(Gd*np.conj(psLp))*dx"),
("dict",
 "    return dict(z=z,t=t,xc=xc,a_w=a_w,rho_b=rho_b,rho_w=rw,D_corr=D_corr,num_direct=numd,\n"
 "                v_formula=vfor,valid=rw>=MASK_REL*rw.max(),p_out=p_out)",
 "    return dict(z=z,t=t,xc=xc,a_w=a_w,rho_b=rho_b,rho_w=rw,D_corr=D_corr,num_direct=numd,\n"
 "                v_formula=vfor,valid=rw>=MASK_REL*rw.max(),p_out=p_out,\n"
 "                rho_b_inc=binsum2(rho_inc*dx),bg_b=binsum2(bg*dx),\n"
 "                V=float(C.V_of(z)),Psc=float(C.Psc_of(z)),SIG=float(C.SIG_of(z)))"),
])

# ============================ stage3_analysis.py ============================
patch("stage3_analysis.py", [
("imports",
 "import numpy as np, pickle, sys\nsys.path.insert(0,'sim')\nd,sig0,w,vz,L=2.5e-6,0.6e-6,2e-6,1.5,0.60",
 "import numpy as np, pickle, sys\nsys.path.insert(0,'sim')\nimport v8_config as C\nW=C.W"),
("mom",
 "mom=((lev*prodb*win).sum(axis=1)+(lev*prodb*(~win)).sum(axis=1)+mp['D_corr'])/w",
 "mom=((lev*prodb*win).sum(axis=1)+(lev*prodb*(~win)).sum(axis=1)+mp['D_corr'])/W"),
("wigprint",
 "{vz*d/zj[k]*1e6:9.2f}",
 "{C.d/C.t_of_z(zj[k])*1e6:9.2f}"),
("wig",
 "wig=res['wiggle_1em2']*0+np.array([vz*d/z for z in zj])   # wiggle scale per plane",
 "wig=np.array([C.d/C.t_of_z(z) for z in zj])   # v8 wiggle scale d/t(z) per plane"),
])

# ============================ stage4_gls.py ============================
patch("stage4_gls.py", [
("imports",
 "from scipy.interpolate import BSpline\nd,sig0,w,vz,L=2.5e-6,0.6e-6,2e-6,1.5,0.60\nGA=4402.0\n"
 "sig_arm=lambda Ng,ga: np.sqrt(Ng)/(np.sqrt(2)*abs(ga))\n"
 "SIG=1.0/np.sqrt(sig_arm(1.684e9,4402.)**-2+sig_arm(4.499e9,4402.)**-2)\nF_SE=0.005",
 "from scipy.interpolate import BSpline\nimport v8_config as C\nd,sig0,L=C.d,C.sig0,C.L; W=C.W"),
("uw",
 "UW=np.where((xf>=-0.5e-3)&(xf<=0.5e-3),1.0,0.0); UW/=UW.sum()\n",
 ""),
("lamf",
 "xc=mp['xc']; t=mp['t']; lam_f=1.224e-3*mp['z']",
 "xc=mp['xc']; t=mp['t']; lam_f=C.lam_f(mp['z'])"),
("momexact",
 "mom_exact=((lev*prodb).sum(axis=1)+mp['D_corr'])/w",
 "mom_exact=((lev*prodb).sum(axis=1)+mp['D_corr'])/W"),
("dict",
 "    return dict(xc=xc,t=t,A=A,B=B,knots=knots,nb=nb,win=win,Mout=Mout,mom_exact=mom_exact,\n"
 "                lo=lo,hi=hi,rho=mp['rho_b'].astype(float),aw=mp['a_w'].astype(float),\n"
 "                Dc=mp['D_corr'],rw=rw,valid=mp['valid'])",
 "    return dict(xc=xc,t=t,A=A,B=B,knots=knots,nb=nb,win=win,Mout=Mout,mom_exact=mom_exact,\n"
 "                lo=lo,hi=hi,rho=mp['rho_b'].astype(float),aw=mp['a_w'].astype(float),\n"
 "                Dc=mp['D_corr'],rw=rw,valid=mp['valid'],\n"
 "                rho_inc=mp['rho_b_inc'].astype(float),bg=mp['bg_b'].astype(float),\n"
 "                V=mp['V'],Psc=mp['Psc'],SIG=mp['SIG'])"),
("momnoisy",
 "def mom_noisy(P,rng,Na):\n"
 "    pmix=(1-F_SE)*P['rho']+F_SE*P['rho'].sum(axis=-1,keepdims=True)*UW[None,:] if P['rho'].ndim==2 else None\n"
 "    rho=P['rho']; pm=(1-F_SE)*rho+F_SE*rho.sum(1,keepdims=True)*UW[None,:]\n"
 "    sigf=np.where(pm>0,(1-F_SE)*rho/np.maximum(pm,1e-300),0.0)\n"
 "    lev=xf[None,:]-P['xc'][:,None]\n"
 "    pv=np.concatenate([pm,np.maximum(1-pm.sum(1,keepdims=True),0)],axis=1)\n"
 "    n=np.stack([rng.multinomial(Na,pv[i]/pv[i].sum())[:NB] for i in range(P['xc'].size)])\n"
 "    ah=sigf*P['aw']+rng.standard_normal(pm.shape)*SIG/np.sqrt(np.maximum(n,1))\n"
 "    dat=(lev*(n/Na)*ah*P['win']).sum(axis=1)\n"
 "    mom=(dat+P['Mout']+P['Dc'])/w\n"
 "    var=((lev**2)*rho*P['win']).sum(axis=1)*SIG**2/(Na*w**2)\n"
 "    return mom,np.sqrt(var)",
 "def mom_noisy(P,rng,Na):\n"
 "    rho=P['rho']\n"
 "    if '_pm' not in P:                                   # v8 mixture, cached per plane\n"
 "        rhoV=P['V']*rho+(1-P['V'])*P['rho_inc'][None,:]  # m12 visibility mixture\n"
 "        base=(1-P['Psc'])*rhoV+P['Psc']*P['bg'][None,:]  # per-plane SE background\n"
 "        sh=lambda s: np.stack([np.interp(xf,xf+s,base[i],left=0,right=0) for i in range(base.shape[0])])\n"
 "        pm=(1-C.F_IMP)*base+0.5*C.F_IMP*(sh(+C.DSCR_IMP)+sh(-C.DSCR_IMP))   # f-impurity split\n"
 "        P['_pm']=pm; P['_sigf']=np.where(pm>0,(1-C.F_IMP)*(1-P['Psc'])*rhoV/np.maximum(pm,1e-300),0.0)\n"
 "    pm=P['_pm']; sigf=P['_sigf']\n"
 "    lev=xf[None,:]-P['xc'][:,None]\n"
 "    pv=np.concatenate([pm,np.maximum(1-pm.sum(1,keepdims=True),0)],axis=1)\n"
 "    n=np.stack([rng.multinomial(Na,pv[i]/pv[i].sum())[:NB] for i in range(P['xc'].size)])\n"
 "    ah=sigf*P['aw']+rng.standard_normal(pm.shape)*P['SIG']/np.sqrt(np.maximum(n,1))\n"
 "    dat=(lev*(n/Na)*ah*P['win']).sum(axis=1)\n"
 "    mom=(dat+P['Mout']+P['Dc'])/W\n"
 "    var=((lev**2)*rho*P['win']).sum(axis=1)*P['SIG']**2/(Na*W**2)\n"
 "    return mom,np.sqrt(var)"),
])

# ============================ stage4_run.py ============================
patch("stage4_run.py", [
("imports",
 "from sim_core import make_psi\nd,sig0,w,vz,L=2.5e-6,0.6e-6,2e-6,1.5,0.60",
 "from sim_core import make_psi\nimport v8_config as C\nd,sig0,L=C.d,C.sig0,C.L; W=C.W"),
("psi",
 "psi,_,_=make_psi(d,sig0,vz)",
 "psi,_,_=make_psi(d,sig0,tmap=C.t_of_z)"),
("traj",
 "    f=lambda z,x: sp(z,np.clip(x,xr[0],xr[-1]),grid=False)/vz",
 "    f=lambda z,x: sp(z,np.clip(x,xr[0],xr[-1]),grid=False)/C.vz_of(z)"),
("momex",
 "        out.append(((lev*prodb).sum(axis=1)+mp['D_corr'])/w)",
 "        out.append(((lev*prodb).sum(axis=1)+mp['D_corr'])/W)"),
("fringe1",
 "({np.mean(tr)/735e-6*100:6.1f}% fringe)",
 "({np.mean(tr)/C.FRINGE*100:6.1f}% fringe)"),
("fringe2",
 "({r/735e-6*100:.3f}% fringe)",
 "({r/C.FRINGE*100:.3f}% fringe)"),
])

# ============================ stage5_alloc.py ============================
patch("stage5_alloc.py", [
("imports",
 "from sim_core import make_psi\nd,sig0,w,vz,L=2.5e-6,0.6e-6,2e-6,1.5,0.60",
 "from sim_core import make_psi\nimport v8_config as C\nd,sig0,L=C.d,C.sig0,C.L; W=C.W"),
("psi",
 "psi,_,_=make_psi(d,sig0,vz)",
 "psi,_,_=make_psi(d,sig0,tmap=C.t_of_z)"),
("traj",
 "    f=lambda z,x: sp(z,np.clip(x,xr[0],xr[-1]),grid=False)/vz",
 "    f=lambda z,x: sp(z,np.clip(x,xr[0],xr[-1]),grid=False)/C.vz_of(z)"),
("sigof",
 "        var=((lev**2)*P['rho']*P['win']).sum(axis=1)*G.SIG**2/(Na*w**2)\n        return np.sqrt(var)",
 "        var=((lev**2)*P['rho']*P['win']).sum(axis=1)*P['SIG']**2/(Na*W**2)\n        return np.sqrt(var)"),
("m16",
 "    G_gain=np.sqrt(20*A.sum())/np.sqrt(A).sum()\n"
 "    print(f\"optimal-allocation gain (this Na, additive model): G = {G_gain:.3f}\")\n"
 "    print(\"optimal shares N_k/N_tot:\", \" \".join(f\"{s:.3f}\" for s in np.sqrt(A)/np.sqrt(A).sum()))",
 "    G_gain=np.sqrt(20*A.sum())/np.sqrt(A).sum()\n"
 "    print(f\"optimal-allocation gain (this Na, additive model): G = {G_gain:.3f}\")\n"
 "    print(\"optimal shares N_k/N_tot:\", \" \".join(f\"{s:.3f}\" for s in np.sqrt(A)/np.sqrt(A).sum()))\n"
 "    # ---- m16: seed-paired (common-random-number) allocation comparison ----\n"
 "    shares_opt=np.sqrt(A)/np.sqrt(A).sum(); shares_uni=np.ones(20)/20.0\n"
 "    NPAIR=12; s0=777000\n"
 "    r_uni=run_alloc(shares_uni,20*Na,NPAIR,s0)\n"
 "    r_opt=run_alloc(shares_opt,20*Na,NPAIR,s0)   # identical seed list: CRN pairing\n"
 "    gain=r_uni/r_opt\n"
 "    store['m16']=dict(Na=Na,shares_opt=shares_opt,r_uni=r_uni,r_opt=r_opt,gain=gain)\n"
 "    pickle.dump(store,open(fn,'wb'))\n"
 "    print(f\"m16 CRN-paired allocation: uniform/optimal RMS ratio = {gain.mean():.3f} +- \"\n"
 "          f\"{gain.std(ddof=1)/np.sqrt(NPAIR):.3f}   ({(gain.mean()-1)*100:+.1f}% benefit)\")"),
])

# ============================ stage6_fieldmap.py ============================
patch("stage6_fieldmap.py", [
("imports",
 "from scipy.ndimage import gaussian_filter1d\nw=2e-6",
 "from scipy.ndimage import gaussian_filter1d\nimport v8_config as C\nW=C.W"),
("lamf",
 "        P=DES[k]; xc=P['xc']; lam_f=1.224e-3*zj[k]",
 "        P=DES[k]; xc=P['xc']; lam_f=C.lam_f(zj[k])"),
("sigof",
 "    return np.sqrt(((lev**2)*P['rho']*P['win']).sum(axis=1)*G.SIG**2/(Na*w**2))",
 "    return np.sqrt(((lev**2)*P['rho']*P['win']).sum(axis=1)*P['SIG']**2/(Na*W**2))"),
])

print("\nALL PATCHES APPLIED")
