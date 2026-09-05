"""Fig 1 (v8.1): vertical photon-curtain apparatus. Script-generated; trajectories and the
screen pattern are the REAL simulation objects (true_traj.pkl, make_psi at z=L), not cartoons.
STATUS NOTE (P8): PARKED alternative to the shipped Gemini schematic — layout label collisions
unresolved; kept in the supplementary for provenance of the computed-screen-profile analysis."""
import numpy as np, pickle, sys
sys.path.insert(0, 'sim')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle, Circle
import v8_config as C
from sim_core import make_psi

BLUE, RED, GOLD, INK = '#2f5fb3', '#b3402f', '#c8a415', '#222222'
plt.rcParams.update({'font.size': 8.5, 'axes.linewidth': 0.8})

TT = pickle.load(open('sim/true_traj.pkl', 'rb'))          # (200, 60) on zeval 5->60 cm
zeval = np.linspace(0.05, 0.60, 60)
psi, _, _ = make_psi(C.d, C.sig0, tmap=C.t_of_z)
xs = np.linspace(-1.6e-3, 1.6e-3, 3000)
rho = np.abs(psi(xs, 0.60))**2; rho /= rho.max()

fig = plt.figure(figsize=(7.2, 5.6))
axA = fig.add_axes([0.06, 0.06, 0.50, 0.90])               # main vertical apparatus
axB = fig.add_axes([0.63, 0.52, 0.345, 0.42])              # curtain zoom
axC = fig.add_axes([0.63, 0.06, 0.345, 0.34])              # screen profile zoom

# ================= (a) main apparatus, z downward =================
a = axA
a.set_xlim(-2.3, 2.6); a.set_ylim(64, -19); a.set_xticks([]); a.set_yticks([])
for sp in a.spines.values(): sp.set_visible(False)
a.add_patch(Circle((0, -13.5), 1.05, fc='#dfe8f5', ec=INK, lw=0.9, zorder=3))
a.text(0, -13.5, r'$^{87}$Rb' + '\nMOT', ha='center', va='center', fontsize=7.6, zorder=4)
a.annotate('', xy=(0, -3.6), xytext=(0, -12.2),
           arrowprops=dict(arrowstyle='-|>', lw=1.0, color=INK))
a.text(0.18, -8.0, 'release,\n11.5 cm fall', fontsize=7.2, va='center')
a.annotate('', xy=(-2.05, 8), xytext=(-2.05, -2),
           arrowprops=dict(arrowstyle='-|>', lw=1.4, color=INK))
a.text(-2.05, 10.8, r'$g$', fontsize=11, ha='center')
for x0, x1 in [(-2.0, -0.09), (0.09, 2.0)]:
    a.add_patch(Rectangle((x0, -0.9), x1-x0, 1.2, fc='0.25', ec='none', zorder=3))
a.text(2.08, -0.3, r'double slit ($d=2.5\,\mu$m,' + '\n' + r'$v_0=1.50$ m/s)', fontsize=7.2, va='center')
for j in range(0, 200, 17):
    a.plot(TT[j]*1e3, zeval*100, lw=0.5, color=BLUE, alpha=0.55, zorder=2)
a.text(-2.25, 26, 'exact Bohmian\ntrajectories\n(transverse axis in mm;\naspect exaggerated)',
       fontsize=6.8, color=BLUE, va='center')
for zc, tag in [(5.0, r'$z_1 = 5$ cm'), (29.6, r'$z_{10}$'), (57.0, r'$z_{20} = 57$ cm')]:
    a.plot([-1.7, 1.7], [zc-0.28, zc-0.28], color=BLUE, lw=1.5, zorder=4)
    a.plot([-1.7, 1.7], [zc+0.28, zc+0.28], color=RED, lw=1.5, zorder=4)
    a.add_patch(Rectangle((1.72, zc-0.75), 0.30, 1.5, fc='#f4e9c8', ec=INK, lw=0.7, zorder=4))
    a.text(2.10, zc, r'$\delta p_\gamma$', fontsize=7.4, va='center')
    a.text(-1.75, zc, tag, fontsize=7.2, ha='right', va='center')
a.annotate('', xy=(1.05, 27.4), xytext=(-1.05, 27.4),
           arrowprops=dict(arrowstyle='<|-|>', lw=0.9, color=INK))
a.text(0, 25.6, r'$x_c$ scanned: 151 positions $\times$ 20 planes' + '\n(one plane active per shot)',
       fontsize=6.9, ha='center')
for zc in np.linspace(0.05, 0.57, 20)*100:
    if abs(zc-5.0) > 1 and abs(zc-29.6) > 1 and abs(zc-57.0) > 1:
        a.plot([-1.45, 1.45], [zc, zc], color='0.82', lw=0.5, ls=(0, (2, 3)), zorder=1)
a.add_patch(Rectangle((-2.0, 59.6), 4.0, 0.55, fc='0.25', ec='none', zorder=3))
a.fill_between(xs*1e3, 60.4, 60.4 + 3.1*rho, color=GOLD, alpha=0.9, zorder=3)
a.text(2.08, 60.0, r'screen ($z=L=60$ cm)', fontsize=7.2, va='center')
a.plot([-1.25, -1.25], [59.3, 64.0], color=INK, lw=0.7, ls=':')
a.plot([ 1.25,  1.25], [59.3, 64.0], color=INK, lw=0.7, ls=':')
a.text(0, 64.6, r'$\pm 1.25$ mm analysis window', fontsize=6.9, ha='center')
a.annotate('', xy=(2.45, 58), xytext=(2.45, 2),
           arrowprops=dict(arrowstyle='-|>', lw=0.8, color='0.4'))
a.text(2.52, 30, r'$v_z(z)$: 1.50 $\to$ 3.74 m/s' + '\n' + r'$T = 0.229$ s', fontsize=6.9,
       rotation=90, va='center', color='0.25')
a.text(-2.25, -17.5, '(a)', fontsize=10, fontweight='bold')

# ================= (b) curtain zoom =================
b = axB
xz = np.linspace(-6, 6, 400)
b.plot(xz, np.exp(-2*xz**2/2.0**2), color=BLUE, lw=1.6, label=r'776.2 nm (blue of D2)')
b.plot(xz, 0.92*np.exp(-2*xz**2/2.0**2), color=RED, lw=1.6, ls='--', label=r'803.5 nm (red of D1)')
b.annotate('', xy=(0.0, 1.55), xytext=(0.0, 2.35),
           arrowprops=dict(arrowstyle='-|>', lw=1.3, color=INK))
b.add_patch(Circle((0.0, 1.38), 0.16, fc=INK, ec='none'))
b.text(0.45, 1.9, 'atom', fontsize=7.4)
b.annotate('', xy=(4.6, 0.42), xytext=(-4.6, 0.12),
           arrowprops=dict(arrowstyle='-|>', lw=1.1, color=GOLD))
b.text(-4.5, -0.13, r'photon in ($\times F{=}100$ passes)', fontsize=6.8, color='#8a7110')
b.text(2.3, 0.55, r'$\delta p_\gamma \propto x_\mathrm{atom}-x_c$', fontsize=7.2, color='#8a7110')
b.axvline(0, color='0.6', lw=0.6, ls=':'); b.text(0.12, -0.32, r'$x_c$', fontsize=8)
b.text(0, 2.75, r'dual-color force balance: $U_{776}+U_{804}=0$;' + '\nwhich-way records add',
       fontsize=6.9, ha='center')
b.set_xlim(-6, 6); b.set_ylim(-0.55, 3.15)
b.set_xticks([]); b.set_yticks([])
b.set_xlabel(r'transverse position (waist $w = 2\,\mu$m)', fontsize=7.2)
b.legend(fontsize=6.3, loc='upper right', frameon=False, borderpad=0.1)
b.text(-5.8, 2.85, '(b)', fontsize=10, fontweight='bold')

# ================= (c) screen pattern zoom =================
c = axC
c.fill_between(xs*1e3, 0, rho, color=GOLD, alpha=0.85)
c.plot(xs*1e3, rho, color='#8a7110', lw=0.8)
pk = xs[np.r_[True, np.diff(np.sign(np.diff(rho))) < 0, True][1:-1] & (rho[1:-1] > 0.25)] if False else None
from scipy.signal import argrelmax
im = argrelmax(rho, order=25)[0]
xc_pk = xs[im][np.argsort(np.abs(xs[im]))][:2]*1e3
x1, x2 = sorted(xc_pk)
c.annotate('', xy=(x2, 1.06), xytext=(x1, 1.06), arrowprops=dict(arrowstyle='<|-|>', lw=0.9, color=INK))
c.text((x1+x2)/2, 1.12, r'$420\,\mu$m fringe', fontsize=7.2, ha='center')
c.axvline(-1.25, color=INK, lw=0.7, ls=':'); c.axvline(1.25, color=INK, lw=0.7, ls=':')
c.set_xlim(-1.6, 1.6); c.set_ylim(0, 1.28)
c.set_yticks([]); c.set_xlabel(r'$x_f$ (mm)', fontsize=7.6)
c.text(-1.52, 1.14, '(c)', fontsize=10, fontweight='bold')
c.text(1.05, 0.92, r'$|\psi(x_f, L)|^2$', fontsize=7.6)

fig.savefig('photon_curtain_fig1_setup.pdf')
fig.savefig('fig1_check.png', dpi=150)
print(f"fig1 saved; fringe between central maxima: {abs(x2-x1)*1e3:.1f} um (anchor 420.3)")
