"""NJP-120031 v8 campaign -- Gate G0: independent verification of decision points D1-D4.

Independent rederivation from CODATA/Steck constants + Table-2 geometry (d, w0, L) only;
no reuse of the plan's arithmetic. Where D3 transforms a v7 chain quantity (Lambda_eff,
kick, budget), the v7 value is taken from the review baseline verify_log.txt (all PASS
there) and the transformation FACTOR is derived here from first principles.
PASS/FAIL vs plan-quoted values (NJP120031_v8_correction_plan.md / v8_plan_numbers.txt).
ASCII output.  2026-07-04.
"""
import numpy as np

# ---- constants (independent) ----
h    = 6.62607015e-34
hbar = 1.054571817e-34
g    = 9.80665
u    = 1.66053906660e-27
M    = 86.909180527 * u          # 87Rb (Steck) -> 1.4431606e-25 kg
Om_e = 7.292115e-5               # Earth sidereal rotation rate (rad/s)
muB  = 9.2740100783e-24
lamD2 = 780.241209686e-9

# ---- Table-2 geometry (unchanged by v8) ----
d, w0, L = 2.5e-6, 2.0e-6, 0.60
lat = 41.3                       # deg (Omaha 41.26 N)
v0  = 1.5                        # D1 default: vz at slit plane

npass = [0, 0]
def ok(name, val, ref, tol=0.02):
    good = abs(val - ref) <= tol * abs(ref)
    npass[0] += good; npass[1] += 1
    print(f"{'PASS' if good else 'FAIL':4s} {name:62s} calc={val:12.6g}  plan={ref:.6g}")
def hdr(s): print("=" * 100 + f"\n{s}\n" + "=" * 100)

hdr("CONSTANTS SANITY")
ok("M(87Rb) (kg)", M, 1.44316e-25, 0.001)
print("NOTE review master_verify.py used M=1.44e-25 (0.22% low); v8 verifier standardizes on Steck value.")

# =============================================================================
hdr("D1 -- ORIENTATION AND LAUNCH (fixes R1): vertical drop, v0(slit) = 1.5 m/s")
hdrop = v0**2 / (2 * g)
ok("drop height above slit plane (cm)", hdrop * 100, 11.5, 0.01)
print(f"NOTE exact drop height = {hdrop*100:.3f} cm; plan's 11.5 is the rounded spec."
      f" Recommend ms quotes v0 = 1.500 m/s primary, h ~ 11.5 cm derived.")
vL = np.sqrt(v0**2 + 2 * g * L)
ok("v(screen) (m/s)", vL, 3.744, 0.001)
T = (vL - v0) / g
ok("total flight T slit->screen (s)", T, 0.2288, 0.001)
fringe = h * T / (M * d)
ok("fringe spacing h T/(M d) (um)", fringe * 1e6, 420, 0.002)
ok("v7 fringe / v8 fringe", 736.2 / (fringe * 1e6), 1.75, 0.005)
v0a = np.sqrt(2 * g * 0.05); vLa = np.sqrt(v0a**2 + 2 * g * L)
Ta = (vLa - v0a) / g; fra = h * Ta / (M * d)
ok("ALT 5-cm drop: v0 (m/s)", v0a, 0.99, 0.001)
ok("ALT 5-cm drop: T (s)", Ta, 0.2631, 0.001)
ok("ALT 5-cm drop: fringe (um)", fra * 1e6, 483, 0.002)

# =============================================================================
hdr("D2 -- CORIOLIS (fixes R2): slit axis x north-south")
# a_E = 2 Om cos(lat) vz(t); integral vz dt over the fall = L  ->  dv_E = 2 Om cos(lat) L,
# independent of the velocity profile.
dvE = 2 * Om_e * np.cos(np.deg2rad(lat)) * L
ok("eastward dv = 2 Om cos(lat) L, lat 41.3 (um/s)", dvE * 1e6, 65.8, 0.005)
ok("worst-case cos(lat) = 1 (um/s)", 2 * Om_e * L * 1e6, 87.5, 0.005)
ok("in-plane residual @ theta = 10 mrad (um/s)", dvE * 1e-2 * 1e6, 0.66, 0.01)
vrec = h / (M * lamD2)
ok("transverse-velocity channel 2 Om v_rec T (um/s)", 2 * Om_e * vrec * T * 1e6, 0.20, 0.03)
print("     derivation check: dv_E = 2 Om cos(lat) * integral(vz dt) = 2 Om cos(lat) * L  [profile-independent]")

# =============================================================================
hdr("D3 -- BEAM CONVENTION (fixes M1) + SIGNS (M2): keep Table 2 (1/e^2), repropagate shapes")
weff = w0 / np.sqrt(2)
ok("w_eff = w0/sqrt(2) (um)", weff * 1e6, 1.41, 0.005)
# D1-coefficient shape factor phi1^2 (d/w)^2 exp(-d^2/2w^2): ratio under w0 -> weff
# (phi1 is per-photon at beam center: peak intensity preserved, so phi1 unchanged)
r_shape = 2 * np.exp(-(d**2 / 2) * (1 / weff**2 - 1 / w0**2))
ok("Lambda shape ratio (w0 -> weff)", r_shape, 0.9157, 0.001)
Leff7 = 0.011321          # review baseline, verify_log PASS
Leff8 = Leff7 * r_shape
ok("Lambda_eff v8", Leff8, 1.04e-2, 0.005)
ok("SNR sqrt(Na Leff) @ Na = 1e4", np.sqrt(1e4 * Leff8), 10.2, 0.005)
ok("Na(SNR = 10) v8", 100 / Leff8, 9647, 0.005)
kick7 = 1.382             # review baseline (per-arm max kick, m/s)
kick8 = kick7 * np.sqrt(2)   # max|d/dx exp(-x^2/weff^2)| = sqrt(2/e)/weff = sqrt(2) x (w0 case)
ok("Delta v_max v8 (m/s)", kick8, 1.95, 0.005)
bud7 = 4.521e-7           # review baseline
ok("two-stage budget v8 (eta*eps)", bud7 / np.sqrt(2), 3.20e-7, 0.005)
# rejected alternative (a): keep declared e^{-x^2/w^2} profile, halve I0 -> phi1/2 -> Lambda/4
ok("ALT(a) SNR @ 1e4", np.sqrt(1e4 * Leff7 / 4), 5.3, 0.01)
ok("ALT(a) Na* multiplier", (100 / (Leff7 / 4)) / (100 / Leff7), 4.0, 1e-9)
print("\n  IMPULSE-CRITERION MARGIN (Jess sign-off item, D3): Dv_max * tau_a(z) vs beam scale")
print("    z[cm]   tau_a[us]   Dv*tau[um]    /w0     /w_eff")
for z in [0.05, 0.077, 0.105, 0.30, 0.57]:
    vz = np.sqrt(v0**2 + 2 * g * z); ta = w0 / vz; dx = kick8 * ta
    print(f"    {z*100:5.1f}   {ta*1e6:8.3f}   {dx*1e6:9.3f}   {dx/w0:5.3f}   {dx/weff:5.3f}")
dx7 = kick7 * (w0 / v0)
print(f"    v7 baseline (uniform tau_a = 1.333 us): Dv*tau = {dx7*1e6:.3f} um -> /w0 = {dx7/w0:.3f}")
print("    -> front plane now EXCEEDS w0 (x1.09); exceeded w_eff at all curtain planes but the last two.")

# =============================================================================
hdr("D4 -- PLANE GRID: per-plane kinematics at the 5 quoted planes")
plan = {0.05:  (1.797, 1.113, 0.0303, 82.4, 75.0),
        0.077: (1.939, 1.031, 0.0448, 55.8, 48.7),
        0.105: (2.076, 0.963, 0.0587, 42.6, 35.7),
        0.30:  (2.852, 0.701, 0.1379, 18.1, 12.5),
        0.57:  (3.665, 0.546, 0.2207, 11.3, 6.6)}
print("    z[cm]   vz[m/s]   tau_a[us]   t(z)[s]    d/t[um/s]   v7 vz0*d/z   signal ratio   row")
for z, (pv, pt, ptt, ps, po) in plan.items():
    vz = np.sqrt(v0**2 + 2 * g * z); ta = w0 / vz * 1e6; tz = (vz - v0) / g
    sig = d / tz * 1e6; old = v0 * d / z * 1e6
    row_ok = (abs(vz - pv) <= 0.005 * pv and abs(ta - pt) <= 0.005 * pt and
              abs(tz - ptt) <= 0.005 * ptt and abs(sig - ps) <= 0.005 * ps and
              abs(old - po) <= 0.005 * po)
    npass[0] += row_ok; npass[1] += 1
    print(f"    {z*100:5.1f}   {vz:7.4f}   {ta:8.4f}   {tz:8.5f}   {sig:9.3f}   {old:9.3f}"
          f"      x{sig/old:5.3f}      {'PASS' if row_ok else 'FAIL'}")
ta5 = w0 / np.sqrt(v0**2 + 2 * g * 0.05); ta57 = w0 / np.sqrt(v0**2 + 2 * g * 0.57)
ok("tau_a ratio first/last plane", ta5 / ta57, 2.04, 0.005)
ok("slit-plane tau_a (us)", w0 / v0 * 1e6, 1.33, 0.005)
print(f"NOTE Table-2 uniform 1.33 us vs true tau_a(z): "
      f"{(ta5*1e6-1.3333)/1.3333*100:+.1f}% at z=5 cm, {(ta57*1e6-1.3333)/1.3333*100:+.1f}% at z=57 cm  [plan: -16% to -59%]")
vz5 = np.sqrt(v0**2 + 2 * g * 0.05)
ok("Lambda(5 cm)/Lambda_table = v0/vz(5 cm)", v0 / vz5, 0.835, 0.002)
ok("5-sigma Na threshold factor (kinematics alone)", vz5 / v0, 1.198, 0.002)
print(f"NOTE opposing-factor heuristic at anchor plane: threshold x{vz5/v0:.3f} (Lambda) vs signal x{(d/((vz5-v0)/g))/(v0*d/0.05):.3f}"
      f" -> naive net x{(vz5/v0)/((d/((vz5-v0)/g))/(v0*d/0.05))**2:.3f} (near-neutral, HEURISTIC ONLY;")
print("      honest recompute through the full matched-filter/tier chain is P2/P7 scope per plan).")
dTdv0 = (v0 / vL - 1) / g
ok("dT/dv0 (ms per m/s)", dTdv0 * 1e3, -61.1, 0.005)
ok("fringe smear (% per cm/s of sigma_v0)", abs(dTdv0) / T * 0.01 * 100, 0.27, 0.02)

# =============================================================================
hdr("P3-PREVIEW NUMBERS QUOTED IN THE PLAN (verify now; full rebuild is P3 scope)")
a1 = muB * 0.5 * 1e-5 / M          # |gF mF| = 1/2, grad B = 1 mG/cm = 1e-5 T/m
ok("mF = +-1 accel @ 1 mG/cm (m/s^2)", a1, 3.21e-4, 0.005)
ok("dv over T (um/s)", a1 * T * 1e6, 73.5, 0.005)
ok("spec {0.1 mG/cm, f = 5%}: ensemble mean (um/s)", a1 * T * 0.1 * 0.05 * 1e6, 0.368, 0.005)
ok("background-gas exposure lo (%)", 0.4 * T / 0.4, 0.23, 0.01)
ok("background-gas exposure hi (%)", 4.0 * T / 0.4, 2.3, 0.01)

print("\n" + "=" * 100)
print(f"SUMMARY: {npass[0]}/{npass[1]} checks PASS")
print("=" * 100)
