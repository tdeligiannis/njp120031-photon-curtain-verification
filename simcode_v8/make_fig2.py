"""Fig 2: what detuning does and does not purchase. Two-level far-detuned scaling, single-pass free-space reference optics.
P8-8b (2026-08-25): panel-(b) dipole kick now uses w_eff = w/sqrt2 as in Eq. (dipole_kick_max) (was the 1/e^2 radius w: curve sqrt2 low, crossover 0.50 -> 0.35 GHz).
All curves from closed forms anchored to frozen constants; no hand-drawn numbers."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

hbar=1.054571817e-34; c=2.998e8; M=1.443e-25
Gam=2*np.pi*6.07e6; Isat=25.03; w=2e-6; Psc=0.0172
weff=w/np.sqrt(2)                # 1/e intensity radius: ALL shape factors, incl. the dipole-kick maximum of Eq. (dipole_kick_max)
om=2*np.pi*(384.230e12)          # omega_L ~ omega_D2 (variation < 2% over range)
vrec=hbar*(2*np.pi/780e-9)/M
sres=3*(780e-9)**2/(2*np.pi)

Dnu=np.logspace(np.log10(0.3e9), np.log10(1e13), 500)   # Delta/2pi
D=2*np.pi*Dnu
phi1=hbar*om*Gam**2/(4*np.pi*w**2*Isat*D)               # Eq. (phi1)
psc1=phi1*Gam/D                                          # p_sc,1 = phi1*Gamma/Delta
ratio=phi1**2/psc1                                       # detuning-independent
kick=np.sqrt(2/np.e)*(hbar/(M*weff))*Psc*D/Gam           # g_a = Psc*Delta/Gamma at fixed Psc; w_eff per Eq. (dipole_kick_max) [P8-8b fix: was w]
stoch=np.sqrt(Psc)*vrec/np.sqrt(3)                       # transverse rms: emission only, isotropic

fig,(a,b)=plt.subplots(1,2,figsize=(10.2,4.1))
# ---- panel (a)
a.loglog(Dnu/1e9, phi1, color='#1f77b4', lw=1.8, label=r'$\varphi_1\propto 1/\Delta$')
a.loglog(Dnu/1e9, psc1, color='#ff7f0e', lw=1.8, label=r'$p_{\mathrm{sc},1}=\varphi_1\Gamma/\Delta\propto 1/\Delta^2$')
a.axhline(ratio[0], color='k', lw=1.8, label=fr'$\varphi_1^2/p_{{\mathrm{{sc}},1}}={ratio[0]:.1e}$ (det.-indep.)')
a.axhline(sres/(2*w**2), color='k', lw=1.4, ls='--', label=fr'ideal bound $\sigma_\mathrm{{res}}/2w^2={sres/(2*w**2):.3f}$')
a.axvline(200, color='0.4', lw=1.0, ls=':')
a.text(200*1.2, 3e-16, r'$\Delta/2\pi=200$ GHz', rotation=90, fontsize=8, color='0.35', va='bottom')
a.axvspan(0.3, 2, color='0.85', zorder=0); a.text(0.75, 3e-16, 'hf structure', fontsize=7, color='0.4', rotation=90, va='bottom')
a.set_xlabel(r'Detuning $\Delta/2\pi$ (GHz)'); a.set_ylabel('dimensionless')
a.set_ylim(1e-16, 0.4); a.set_xlim(0.3, 1e4)
a.legend(fontsize=7.5, loc='upper right'); a.set_title('(a) information per scattering is detuning-independent', fontsize=9.5)
# ---- panel (b)
b.loglog(Dnu/1e9, kick, color='#2ca02c', lw=1.8,
         label=r'deterministic dipole kick $\sqrt{2/e}\,(\hbar/Mw_\mathrm{eff})\,P_\mathrm{sc}\Delta/\Gamma$')
b.axhline(stoch, color='#9467bd', lw=1.8, label=fr'rms transverse recoil $\sqrt{{P_\mathrm{{sc}}/3}}\,\hbar k/M$')
xcross=np.interp(stoch, kick, Dnu)/1e9
b.plot(xcross, stoch, 'ko', ms=5); b.annotate(fr'crossover $\approx{xcross:.2f}$ GHz', (xcross, stoch),
        textcoords='offset points', xytext=(10,-16), fontsize=8)
b.axvline(200, color='0.4', lw=1.0, ls=':')
b.axvspan(0.3, 2, color='0.85', zorder=0)
b.set_xlabel(r'Detuning $\Delta/2\pi$ (GHz)'); b.set_ylabel('transverse velocity (m/s)')
b.set_xlim(0.3, 1e4); b.set_ylim(1e-5, 3)
b.legend(fontsize=7.5, loc='upper left'); b.set_title(r'(b) at fixed $P_\mathrm{sc}=0.017$: character, not information', fontsize=9.5)
fig.tight_layout()
fig.savefig('photon_curtain_fig2_dispersive.pdf'); fig.savefig('fig2_check.png', dpi=140)
print(f"phi1(200 GHz)={np.interp(200e9,Dnu,phi1):.3e}  psc1(200 GHz)={np.interp(200e9,Dnu,psc1):.3e}")
print(f"ratio={ratio[0]:.3e}  ideal={sres/(2*w**2):.4f}  kick(200 GHz)={np.interp(200e9,Dnu,kick):.3f} m/s  stoch={stoch*1e3:.3f} mm/s  crossover={xcross:.2f} GHz")
