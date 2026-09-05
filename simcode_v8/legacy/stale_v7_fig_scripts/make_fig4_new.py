"""Fig 4: forward-modeled reconstruction at Config E through the GLS pipeline (Tier-2 demo at Na=1e7)."""
import numpy as np, pickle, sys
sys.path.insert(0,'sim')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import stage4_gls as G
plt.rcParams.update({'font.size':9,'axes.labelsize':10,'axes.titlesize':10,'figure.dpi':150})
F=pickle.load(open('sim/fig4_data.pkl','rb'))
TT=pickle.load(open('sim/true_traj.pkl','rb'))
DES=pickle.load(open('sim/stage4_design.pkl','rb'))
zeval=np.linspace(0.05,0.60,60); zj=np.linspace(0.05,0.57,20)
Vt0=G.gls_fit(DES[0],DES[0]['mom_exact']-DES[0]['rw'],np.ones(DES[0]['xc'].size))
# field bands at plane 0
bands={}
for Na in [1e6,1e7]:
    vs=[]
    for r in range(12):
        rng=np.random.default_rng(90000+r)
        mom,sg=G.mom_noisy(DES[0],rng,int(Na))
        vs.append(G.gls_fit(DES[0],mom-DES[0]['rw'],sg))
    bands[Na]=np.array(vs)
fig,ax=plt.subplots(2,2,figsize=(9.5,7.2))
sub=slice(0,200,3)
for a,(T,tt) in zip(ax[0],[ (TT,'(a) true Bohmian trajectories'), (F['rep_1e7'][1],'(b) reconstructed, $N_a=10^7$/point (GLS)') ]):
    a.plot(zeval*100,np.array(T)[sub].T*1e6,lw=0.35,color='C0',alpha=0.7)
    a.set_xlabel('z (cm)'); a.set_ylabel(r'x ($\mu$m)'); a.set_title(tt,loc='left'); a.set_ylim(-650,650)
for zc in zj: ax[0,1].axvline(zc*100,color='g',ls=(0,(1,3)),lw=0.5,alpha=0.5)
P0=DES[0]; core=P0['rw']>=0.05*P0['rw'].max()
a=ax[1,0]
a.plot(P0['xc'][core]*1e6,Vt0[core]*1e3,'k-',lw=1.4,label='true (formula) field')
for Na,c in [(1e6,'C1'),(1e7,'C2')]:
    m=bands[Na].mean(0); s=bands[Na].std(0)
    a.fill_between(P0['xc'][core]*1e6,(m-s)[core]*1e3,(m+s)[core]*1e3,color=c,alpha=0.25)
    a.plot(P0['xc'][core]*1e6,m[core]*1e3,c,lw=0.9,label=f'GLS, $N_a=10^{{{int(np.log10(Na))}}}$ ($\\pm1\\sigma$)')
a.set_xlabel(r'$x_c$ ($\mu$m)'); a.set_ylabel(r'$v_x$ (mm/s)'); a.set_title('(c) velocity field, z = 5.0 cm',loc='left'); a.legend(fontsize=7.5,frameon=False)
a=ax[1,1]
for Na,c in [(1e6,'C1'),(1e7,'C2')]:
    cur=F[Na]['rms_curves']; m=cur.mean(0); se=cur.std(0)/np.sqrt(20)
    a.fill_between(zeval*100,(m-se)*1e6,(m+se)*1e6,color=c,alpha=0.25)
    a.plot(zeval*100,m*1e6,c,lw=1.2,label=f'$N_a=10^{{{int(np.log10(Na))}}}$/point')
a.plot(zeval*100,F['clean']['rms']*1e6,'k--',lw=1,label='noiseless-pipeline floor')
a.axhline(0.33*735,color='gray',ls=':',lw=1); a.text(6,0.33*735*1.08,'33% of fringe',fontsize=7.5,color='gray')
a.set_xlabel('z (cm)'); a.set_ylabel(r'trajectory RMS ($\mu$m)'); a.set_title('(d) RMS deviation vs propagation',loc='left'); a.set_yscale('log'); a.legend(fontsize=7.5,frameon=False)
plt.tight_layout(); plt.savefig('photon_curtain_fig4_reconstruction.pdf'); plt.savefig('photon_curtain_fig4_reconstruction.png',dpi=200)
print("fig4 saved")
