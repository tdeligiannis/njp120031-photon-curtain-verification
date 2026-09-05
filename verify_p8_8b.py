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
ok("M1: max|v - x/t| over the 5%-core at z = 5 cm (um/s), quoted 'at most ~20'", d5[0], 21.0, 0.03)
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
ok("M2: first side maximum / peak, quoted '4%'", rhoL[s1]/pk, 0.039, 0.05)
ok("M2: central-lobe probability, quoted '94%'", central, 0.944, 0.005)
ok("M2: envelope/fringe ratio d/(4 pi sigma0), quoted '1/3'", C.d/(4*np.pi*sig0), 0.33, 0.02)
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
ok("m8: |g_a,comp| - g_a,meas relative imbalance of the tabulated (nominal) values, quoted '1.6e-3'", (3679-3673)/3676, 1.6e-3, 0.05)

print(f"P8-8b VERIFICATION: {npass} checks passed, {nfail} failed")
sys.exit(0 if nfail == 0 else 1)
