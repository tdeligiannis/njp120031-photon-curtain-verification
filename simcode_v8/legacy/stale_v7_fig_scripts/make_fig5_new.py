"""Fig 5: the staged measurement program — trajectory RMS vs Na and field-map accuracy vs Na."""
import numpy as np, pickle, sys
sys.path.insert(0,'sim')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9,'axes.labelsize':10,'figure.dpi':150})
F=pickle.load(open('sim/fig4_data.pkl','rb'))
S6=pickle.load(open('sim/stage6_results.pkl','rb'))
R4=pickle.load(open('sim/stage4_results.pkl','rb')); snr=R4['snr_1e4']
fig,ax=plt.subplots(1,2,figsize=(9.5,3.9))
a=ax[0]
Nas=[1e4,1e5,1e6,1e7]
m=[F[Na]['rms_curves'][:,-1].mean() for Na in Nas]; se=[F[Na]['rms_curves'][:,-1].std()/np.sqrt(20) for Na in Nas]
a.errorbar(Nas,np.array(m)*1e6,yerr=np.array(se)*1e6,fmt='o-',color='C2',lw=1.2,ms=4,capsize=2,label='GLS pipeline (this work)')
raw=[(1e4,22.114e-3),(1e5,9.842e-3),(1e6,5.534e-3)]
a.plot([r[0] for r in raw],[r[1]*1e6 for r in raw],'s--',color='C3',ms=4,lw=1,label='direct moment pipeline')
a.plot(Nas,np.array(m[0])*1e6*np.sqrt(1e4/np.array(Nas)),color='gray',ls=':',lw=0.8)
a.text(2.5e6,m[0]*1e6*np.sqrt(1e4/2.5e6)*1.25,r'$N_a^{-1/2}$',fontsize=8,color='gray')
a.axhline(0.33*735,color='k',ls='-.',lw=0.8); a.text(1.3e4,0.33*735*1.15,'33% of fringe',fontsize=7.5)
a.axhline(4.4,color='k',ls='--',lw=0.8); a.text(1.3e4,5.4,'noiseless-pipeline floor',fontsize=7.5)
a.set_xscale('log'); a.set_yscale('log'); a.set_xlabel(r'$N_a$ per grid point'); a.set_ylabel(r'trajectory RMS at $z=L$ ($\mu$m)')
a.set_title('(a) Tier 2: trajectory reconstruction',loc='left'); a.legend(fontsize=7.5,frameon=False)
a=ax[1]
cols=['C0','C1','C4']; zlbl=['5.0','7.7','10.5']
for k in [0,1,2]:
    Ns=[3e5,1e6,3e6,1e7,3e7]
    mm=[np.mean(S6[(k,Na)]['rel'])*100 for Na in Ns]; ss=[np.std(S6[(k,Na)]['rel'])/np.sqrt(12)*100 for Na in Ns]
    a.errorbar(Ns,mm,yerr=ss,fmt='o-',color=cols[k],ms=3.5,lw=1,capsize=2,label=f'z = {zlbl[k]} cm')
    a.axvline(1e4*(5/snr[k])**2,color=cols[k],ls=':',lw=0.8)
a.axhline(20,color='gray',ls='--',lw=0.8); a.text(3.4e5,21.5,'20%',fontsize=7.5,color='gray')
a.set_xscale('log'); a.set_yscale('log'); a.set_xlabel(r'$N_a$ per grid point'); a.set_ylabel('field relative error (%)')
a.set_title(r'(b) Tier 1: field maps (dotted: 5$\sigma$ detection)',loc='left'); a.legend(fontsize=7.5,frameon=False)
plt.tight_layout(); plt.savefig('photon_curtain_fig5_program.pdf'); plt.savefig('photon_curtain_fig5_program.png',dpi=200)
print("fig5 saved")
