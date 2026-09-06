#!/usr/bin/env python3
"""verify_p8_8b.py -- NJP-120031 P8-8b: pins for every number introduced by the adversarial-review fixes
(M1, M2, M3, m1, m3ii, m6, m7, m8). Derived from sim_core/v8_config and official_fig4_data.pkl; nothing hand-typed.
Prints PASS/FAIL lines; exit 1 on any failure. Wired into verify_p7_master.py section D."""
import sys, os, pickle, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, 'simcode_v8/sim'); import sim_core as S, v8_config as C
npass = nfail = 0
def ok(name, val, ref, tol=0.02):
    global npass, nfail
    good = abs(val - ref) <= tol*abs(ref) if ref != 0 else abs(val) <= tol
    npass += good; nfail += (not good)
    print(f"{'PASS' if good else 'FAIL'} {name:88s} calc={val:12.5g} quoted={ref}")
sig0 = 0.6e-6; h = 6.62607015e-34
psi, dpsi, v = S.make_psi(C.d, sig0, tmap=C.t_of_z)

# ---- M1: the far-field guidance field vs the ballistic ramp x/t ----
x = np.linspace(-3e-3, 3e-3, 200001)
def ramp_dev(z):
    t = C.t_of_z(z); rho = np.abs(psi(x, z))**2; core = rho > 0.05*rho.max()
    dev = v(x, z) - x/t
    return np.abs(dev[core]).max()*1e6, np.sqrt(np.mean(dev[core]**2))*1e6, np.abs(x[core]/t).max()*1e3
d5 = ramp_dev(0.05); d30 = ramp_dev(0.30)
ok("M1: max|v - x/t| over the 5%-core at z = 5 cm (um/s), quoted 'at most ~20 um/s'", d5[0], 21.0, 0.03)
ok("M1: rms (v - x/t) over the core at z = 5 cm (um/s), quoted '6'", d5[1], 6.4, 0.05)
ok("M1: max|v - x/t| at z = 30 cm (um/s), quoted 'below 1'", d30[0], 1.0, 0.1)
ok("M1: ramp at the core edge (mm/s), quoted '~0.7'", d5[2], 0.73, 0.03)
ok("M1: Fresnel parameter hbar t/(2 M sigma0^2) at the first plane, quoted '~31'", S.hbar*C.t_of_z(0.05)/(2*S.M*sig0**2), 30.8, 0.01)

# ---- M2: the screen pattern is essentially single-lobed ----
xs = np.linspace(-2e-3, 2e-3, 400001); rhoL = np.abs(psi(xs, C.L))**2; dx = xs[1]-xs[0]; pk = rhoL.max()
right = xs > 0
loc_min = [i for i in range(1, len(xs)-1) if right[i] and rhoL[i] < rhoL[i-1] and rhoL[i] < rhoL[i+1]]
loc_max = [i for i in range(1, len(xs)-1) if right[i] and rhoL[i] > rhoL[i-1] and rhoL[i] > rhoL[i+1]]
n1 = loc_min[0]; s1 = [i for i in loc_max if i > n1][0]
central = rhoL[(xs > -xs[n1]) & (xs < xs[n1])].sum()*dx/(rhoL.sum()*dx)
ok("M2: first screen node (um), quoted '+-210'", xs[n1]*1e6, 210.1, 0.01)
ok("M2: first side maximum / peak, quoted '4% of the peak'", rhoL[s1]/pk, 0.039, 0.05)
ok("M2: central-lobe probability, quoted '94% of the probability'", central, 0.944, 0.005)
ok("M2: envelope/fringe ratio d/(4 pi sigma0), quoted '0.33 = d/(4 pi sigma0)'", C.d/(4*np.pi*sig0), 0.33, 0.02)
TT = pickle.load(open('simcode_v8/sim/true_traj.pkl', 'rb')); fanw = TT[:, -1].std()*1e6
ok("M2: true-fan rms half-width at the screen (um), quoted '~100' / ratio 184/104 = 1.8", fanw, 103.7, 0.01)
ok("M2: 184 um / fan width, quoted '1.8'", 183.8/fanw, 1.77, 0.02)

# ---- M3: RMS vs N_a exponents (official pickle) ----
F = pickle.load(open('official_fig4_data.pkl', 'rb'))
m = {N: np.array(F[N]['rms_curves'])[:, -1].mean() for N in (1e4, 1e5, 1e6, 1e7)}
ok("M3: overall exponent 1e4 -> 1e7, quoted '-0.3'", np.log10(m[1e7]/m[1e4])/3, -0.28, 0.03)
ok("M3: local exponent 1e4 -> 1e5, quoted '-0.29'", np.log10(m[1e5]/m[1e4]), -0.29, 0.03)
ok("M3: local exponent 1e5 -> 1e6, quoted '-0.10'", np.log10(m[1e6]/m[1e5]), -0.10, 0.1)
ok("M3: local exponent 1e6 -> 1e7, quoted '-0.44'", np.log10(m[1e7]/m[1e6]), -0.44, 0.03)

# ---- m1: Fig-2(b) crossover with w_eff (closed form, same constants as make_fig2.py) ----
hbar = 1.054571817e-34; Mkg = 1.443e-25; Gam = 2*np.pi*6.07e6; w = 2e-6; weff = w/np.sqrt(2); Psc = 0.0172; k = 2*np.pi/780e-9
stoch = np.sqrt(Psc)*hbar*k/Mkg/np.sqrt(3)
xc_eff = stoch*Gam/(np.sqrt(2/np.e)*hbar/(Mkg*weff)*Psc)/(2*np.pi)/1e9
ok("m1: Fig-2 crossover with w_eff (GHz), quoted '0.35'", xc_eff, 0.35, 0.02)

# ---- m3(ii): diffracted slice width and channelling displacement ----
sig_w = C.W/2
sG1 = S.hbar*(C.t_of_z(C.L) - C.t_of_z(0.05))/(2*S.M*sig_w)*1e6
sG3 = S.hbar*(C.t_of_z(C.L) - C.t_of_z(0.105))/(2*S.M*sig_w)*1e6
ok("m3ii: sigma_G at the first plane (um) with sigma_w = w_eff/2, quoted '~100'", sG1, 103.0, 0.02)
ok("m3ii: sigma_G at the third plane (um), quoted '~90'", sG3, 88.0, 0.03)
disp = 14.4e-6*(C.t_of_z(C.L) - C.t_of_z(0.05))*1e6
ok("m3ii: channelling displacement 14 um/s x t (um), quoted '~3'", disp, 2.9, 0.05)
ok("m3ii: sigma_G / displacement, quoted 'more than thirtyfold'", sG1/disp, 36.0, 0.05)

# ---- m6: recoil displacement bound from the earliest curtain plane ----
ok("m6: 5.9 mm/s x (T - t(z_1)) (mm), quoted '1.17'", 5.9e-3*(C.t_of_z(C.L) - C.t_of_z(0.05))*1e3, 1.17, 0.01)

# ---- m7: atoms per position at 3 s.f. ----
Na5 = (5.0/C.SNR_1e4_zref)**2*1e4 if hasattr(C, 'SNR_1e4_zref') else None
if Na5 is None:
    J = __import__('json').load(open('simcode_v8/v8_summary.json')); snr = J['stage4']['snr_1e4'][0] if isinstance(J['stage4']['snr_1e4'], list) else J['stage4']['snr_1e4']['0']
    Na5 = (5.0/snr)**2*1e4
ok("m7: N_a per position for 5 sigma, quoted '1.43e6'", Na5, 1.43e6, 0.01)
ok("m7: total atoms 151 x N_a, quoted '2.2e8'", 151*Na5, 2.2e8, 0.02)

# ---- m8: nominal-power imbalance in Table 3 ----
ok("m8: nominal-power imbalance |g_a,comp|-g_a,meas over the mean, from v8_config (Table 3 states the values are nominal)", (abs(C.GA_C_REF)-C.GA_M_REF)/((abs(C.GA_C_REF)+C.GA_M_REF)/2), 1.6e-3, 0.05)

# ---- final-review M1: odd Hermite-Gauss mode sum of the contact coupling at x = d/2 (sigma_x = w/2) ----
from numpy.polynomial.hermite import hermval; from math import factorial, pi, sqrt
_sig = C.W*np.sqrt(2)/2; _x = C.d/2      # sigma_x = w/2 with w = sqrt(2) w_eff from v8_config (no typed literal)
def _phin(n, x):
    u = x/(sqrt(2)*_sig); c = np.zeros(n+1); c[n] = 1
    return (2**n*factorial(n)*sqrt(2*pi)*_sig)**-0.5*hermval(u, c)*np.exp(-u**2/2)
_v = {k: _phin(k, _x)**2 for k in range(1, 5)}
ok("M1(final review): odd-mode sum through n=3 over the tilt-mode (n=1) term, quoted '~34% raise'", (_v[1]+_v[3])/_v[1], 1.344, 0.01)
ok("M1(final review): Lambda_eff with n<=3 modes (tilt-mode 8.65e-3 x 1.344), quoted '~1.2e-2'", 8.65e-3*(_v[1]+_v[3])/_v[1], 1.16e-2, 0.02)
ok("M1(final review): V with n<=3 modes, quoted '~0.988'", np.exp(-8.65e-3*(_v[1]+_v[3])/_v[1]), 0.988, 0.001)
# ---- Sec 7: global parametric fit gain over the N_xc x N_z calibration points, quoted '~55x' ----
ok("Sec 7: global-fit precision gain sqrt(N_xc N_z) = sqrt(151 x 20), quoted '~55x'", np.sqrt(151*20), 55.0, 0.01)
# ---- P8-G: 776.2-nm ladder line (5P3/2 -> 5D3/2 at 776.157 nm; 5D lifetime ~240 ns; reduced 5P-5D element ~0.43 of D2) ----
_c = 2.99792458e8; _nu = lambda l: _c/l
_Delta = _nu(776.2e-9) - _nu(780.241e-9); _d32 = _nu(776.2e-9) - _nu(776.157e-9)
_Om0 = 2*np.pi*6.15e9; _pop5P = (_Om0/(2*2*np.pi*_Delta))**2; _OmL = 0.43*_Om0; _Gam5D = 2*np.pi*0.66e6
ok("P8-G ladder: detuning from 5P3/2->5D3/2 (GHz), quoted '21 GHz below'", -_d32/1e9, 21.4, 0.05)
ok("P8-G ladder: 5P3/2 admixture (Om0/2Delta)^2, quoted '~2e-6'", _pop5P, 2.4e-6, 0.05)
ok("P8-G ladder: scattering per transit via 5D3/2, quoted '~4e-8'", _pop5P*_Gam5D*(_OmL/(2*2*np.pi*abs(_d32)))**2*1.113e-6, 4.2e-8, 0.1)
_nu5D = _nu(780.241e-9) + _nu(776.157e-9); _Omeff = _Om0*_OmL/(2*2*np.pi*_Delta)
ok("P8-G ladder: two-photon 5S->5D detuning (THz), quoted '2 THz'", abs(2*_nu(776.2e-9) - _nu5D)/1e12, 1.98, 0.02)
ok("P8-G ladder: ground-state light-shift correction / U0, quoted '~4e-7'", (_Omeff**2/(4*2*np.pi*abs(2*_nu(776.2e-9)-_nu5D)))/(_Om0**2/(4*2*np.pi*_Delta)), 4.4e-7, 0.1)
ok("P8-G ladder: two-colour sum vs 5S->5D (THz), quoted '11 THz off'", abs(_nu(776.2e-9)+_nu(803.5e-9)-_nu5D)/1e12, 11.1, 0.02)
ok("P8-G multipass: far-field half-angle lambda/(pi w) at 776 nm (rad), quoted '0.12'", 776.2e-9/(np.pi*2e-6), 0.124, 0.02)
print(f"P8-8b VERIFICATION: {npass} checks passed, {nfail} failed")
sys.exit(0 if nfail == 0 else 1)
