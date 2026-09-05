"""Fig 3: (F, N_a) regime map at fixed P_sc = 0.017. All numbers derived from frozen Config E constants."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba

# Frozen constants (aprime_constants_frozen.txt)
Lam100   = 8.65e-3           # Lambda_eff at F=100, Config E (v8 corrected couplings)
Lam_sp   = 8.6e-5            # single-pass free-space reference (v8; Table 2 ref block)
Na_star, F_star = 1e4, 100   # Configuration E operating point
SNR_star = np.sqrt(Na_star*Lam100)

Lam = lambda F: Lam100*F/100.0
Na_floor = lambda F: 100.0/Lam(F)          # SNR = 10 statistics floor
F_coh = 100.0/Lam100                        # Lambda_eff = 1 boundary

F  = np.logspace(0, np.log10(3e4), 400)
Na = np.logspace(2, 8, 400)
FF, NN = np.meshgrid(F, Na)
SNR = np.sqrt(NN*Lam(FF))

fig, ax = plt.subplots(figsize=(7.2, 5.0))
# regions
viable = (SNR >= 10) & (Lam(FF) < 1)
weak_stats = SNR < 10
strong = Lam(FF) >= 1
ax.contourf(FF, NN, viable.astype(float), levels=[0.5,1.5], colors=[to_rgba('#2ca02c',0.18)])
ax.contourf(FF, NN, weak_stats.astype(float), levels=[0.5,1.5], colors=[to_rgba('0.5',0.25)])
ax.contourf(FF, NN, strong.astype(float), levels=[0.5,1.5], colors=[to_rgba('#d62728',0.20)])
# SNR contours
cs = ax.contour(FF, NN, SNR, levels=[1,3,10,30,100], colors='k', linewidths=[0.7,0.7,1.6,0.7,0.7])
ax.clabel(cs, fmt=lambda v: f"SNR = {v:g}", fontsize=8, inline=True)
# coherence boundary
ax.axvline(F_coh, color='#d62728', lw=1.6)
ax.text(F_coh*1.12, 3e3, r'$\Lambda_\mathrm{eff}=1$', color='#d62728', rotation=90, va='bottom', fontsize=9)
# markers
ax.plot(F_star, Na_star, marker='*', ms=17, mfc='#ffdf00', mec='k', mew=0.9, ls='none', zorder=6,
        label=fr'Configuration E ($F=100$, $N_a=10^4$, SNR $\approx {SNR_star:.1f}$)')
Na30 = Na_floor(30.0)
ax.plot(30, Na30, marker='o', ms=8, mfc='none', mec='k', mew=1.2, ls='none', zorder=6,
        label=fr'$F=30$ variant at SNR $=10$ ($N_a\approx {Na30/1e4:.1f}\times10^4$)')
Na_sp10 = 100.0/Lam_sp
ax.plot(1, Na_sp10, marker='s', ms=7, mfc='w', mec='k', mew=1.2, ls='none', zorder=6,
        label=fr'single pass, free-space optics ($N_a\approx1.2\times10^6$)')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(1, 3e4); ax.set_ylim(1e2, 1e8)
ax.set_xlabel(r'Multipass number $F$   (detunings rescaled $\propto\sqrt{F}$, $P_\mathrm{sc}=0.017$ fixed)')
ax.set_ylabel(r'Atoms per grid point $N_a$')
# secondary top axis: Lambda_eff
top = ax.secondary_xaxis('top', functions=(lambda f: Lam(f), lambda l: 100*l/Lam100))
top.set_xlabel(r'$\Lambda_\mathrm{eff}(F)$')
ax.text(2.2, 2.5e2, 'insufficient statistics', color='0.35', fontsize=9)
ax.text(1.6e3, 3e6, 'viable', color='#1a7a1a', fontsize=10)
ax.text(1.35e4, 3e5, 'strong\nmeasurement', color='#a02020', fontsize=9, ha='center')
ax.legend(loc='upper right', fontsize=8, framealpha=0.95)
fig.tight_layout()
fig.savefig('photon_curtain_fig3_regime.pdf'); fig.savefig('fig3_check.png', dpi=140)
# print the derived anchor numbers for the record
print(f"F_coh(Lam=1) = {F_coh:.0f}")
print(f"Na_floor(F=1) = {Na_floor(1):.3e} ; ref-block single-pass Na@10 = {Na_sp10:.3e}")
print(f"Na_floor(F=30) = {Na30:.3e} ; Na_floor(F=100) = {Na_floor(100):.3e} ; SNR(star) = {SNR_star:.2f}")
