"""NJP-120031 v8 -- P1 (Gate G1): kinematics module verification.
Reproduces EVERY kinematic number entering the Sec.-2 module, Table 2 (kinematic rows),
and the kinematics-driven statements swept in P1. Grid planes are the sim's actual
zj = linspace(0.05, 0.57, 20); the plan's illustrative 30-cm row is realized as the
nearest actual plane (29.63 cm) -- logged as micro-deviation P1-a. ASCII output."""
import numpy as np

h = 6.62607015e-34; hbar = 1.054571817e-34; kB = 1.380649e-23
g = 9.80665; u = 1.66053906660e-27
M = 86.909180527 * u
lamD2 = 780.241209686e-9
d, w0, L = 2.5e-6, 2.0e-6, 0.60
v0 = 1.5

npass = [0, 0]
def ok(name, val, ref, tol=0.005):
    good = abs(val - ref) <= tol * abs(ref)
    npass[0] += good; npass[1] += 1
    print(f"{'PASS' if good else 'FAIL':4s} {name:64s} calc={val:12.6g}  ms={ref:.6g}")
def hdr(s): print("=" * 104 + f"\n{s}\n" + "=" * 104)

hdr("SEC. 2 MODULE / TABLE-2 KINEMATIC ROWS")
ok("drop height v0^2/2g (cm)  [ms: '~11.5']", v0**2/(2*g)*100, 11.5, 0.01)
vL = np.sqrt(v0**2 + 2*g*L)
ok("v_z(L) (m/s)", vL, 3.744, 0.001)
T = (vL - v0)/g
ok("T = t(L) (s)", T, 0.2288, 0.001)
ok("fringe hT/(Md) (um)", h*T/(M*d)*1e6, 420, 0.002)
ok("lambda_dB at slit (nm)", h/(M*v0)*1e9, 3.06, 0.002)
ok("kinetic 'temperature' at slit (mK) [ms '~12']", 0.5*M*v0**2/kB*1e3, 12, 0.03)
ok("coherence velocity hbar/(Md) (mm/s)  [row unchanged]", hbar/(M*d)*1e3, 0.29, 0.02)
ok("recoil velocity (mm/s)  [row unchanged]", h/(M*lamD2)*1e3, 5.9, 0.005)

hdr("PER-PLANE TABLE (tab:kinematics) -- ACTUAL GRID PLANES k = 0,1,2,9,19 of linspace(0.05,0.57,20)")
zg = np.linspace(0.05, 0.57, 20)
ms_rows = {0: (5.00, 1.797, 1.113, 0.0303, 82.4),
           1: (7.74, 1.941, 1.030, 0.0450, 55.6),
           2: (10.47, 2.075, 0.964, 0.0586, 42.7),
           9: (29.63, 2.839, 0.704, 0.1366, 18.3),
           19: (57.00, 3.665, 0.546, 0.2207, 11.3)}
print("    k   z[cm]    vz[m/s]   tau_a[us]   t(z)[s]    d/t[um/s]   row")
for k, (mz, mv, mt, mtt, ms_) in ms_rows.items():
    z = zg[k]; vz = np.sqrt(v0**2 + 2*g*z); ta = w0/vz*1e6; tz = (vz - v0)/g; sig = d/tz*1e6
    row_ok = (abs(z*100 - mz) <= 0.005*mz and abs(vz - mv) <= 0.001*mv and abs(ta - mt) <= 0.002*mt
              and abs(tz - mtt) <= 0.005*mtt and abs(sig - ms_) <= 0.005*ms_)
    npass[0] += row_ok; npass[1] += 1
    print(f"   {k:2d}  {z*100:6.2f}   {vz:7.4f}   {ta:8.4f}   {tz:8.5f}   {sig:9.3f}    {'PASS' if row_ok else 'FAIL'}")
ta_first = w0/np.sqrt(v0**2 + 2*g*zg[0])*1e6; ta_last = w0/np.sqrt(v0**2 + 2*g*zg[19])*1e6
ok("Sec.-3 statement: tau_a at z_ref = 5 cm (us)", ta_first, 1.11, 0.005)
ok("Sec.-3 statement: tau_a at last plane (us)  [ms '0.55']", ta_last, 0.55, 0.01)
print("NOTE plan's illustrative 30-cm row -> nearest actual plane 29.63 cm: vz 2.852->2.839, tau 0.701->0.704,")
print("     t 0.1379->0.1366, d/t 18.1->18.3.  [micro-deviation P1-a, logged]")

hdr("LONGITUDINAL-SPREAD BUDGET + SPEC")
dTdv0 = (v0/vL - 1)/g
ok("dT/dv0 (ms per m/s)", dTdv0*1e3, -61.1, 0.005)
ok("fringe smear (% per cm/s of sigma_v0)", abs(dTdv0)/T*0.01*100, 0.27, 0.02)
ok("smear at spec sigma_v0 = 1 cm/s (%)  [ms 'below 0.3%']", abs(dTdv0)/T*0.01*100, 0.27, 0.02)

hdr("KINEMATICS-DRIVEN STATEMENT SWEEP (Secs. 6-9)")
t_fall = v0/g
ok("fall time release -> slits (s)  [ms '0.153']", t_fall, 0.153, 0.005)
ok("irreducible free-fall piece (ms)  [ms '382']", (t_fall + T)*1e3, 382, 0.005)
cyc_lo = 0.5 + 0.008 + 0.020 + (t_fall + T) + 0.025 + 0.050
cyc_hi = 2.0 + 0.012 + 0.080 + (t_fall + T) + 0.050 + 0.200
print(f"NOTE cycle-time window: {cyc_lo:.2f}--{cyc_hi:.2f} s -> footnote's '~1--2.7 s' statement SURVIVES unchanged.")
occ_lo = 1e2*ta_last*1e-6/10e-3; occ_hi = 1e3*ta_first*1e-6/10e-3
ok("multi-atom occupancy low  [ms rounded '0.005']", occ_lo, 0.005, 0.10)
ok("multi-atom occupancy high [ms rounded '0.1']", occ_hi, 0.1, 0.12)
ok("recoil displacement from the earliest curtain plane 5.9 mm/s x (T - t(0.05)) (mm)  [ms '1.17'; P8-8b m6]", h/(M*lamD2)*(T - (np.sqrt(v0**2 + 2*g*0.05) - v0)/g)*1e3, 1.17, 0.005)
sig_first = d/((np.sqrt(v0**2+2*g*zg[0])-v0)/g)*1e6
sig_last = d/((np.sqrt(v0**2+2*g*zg[19])-v0)/g)*1e6
ok("signal scale first plane (um/s)  [for P2 Sec.-7 rewrite]", sig_first, 82.4, 0.005)
ok("signal scale last plane (um/s)   [for P2 Sec.-7 rewrite]", sig_last, 11.3, 0.005)

print("\n" + "=" * 104)
print(f"SUMMARY: {npass[0]}/{npass[1]} checks PASS")
print("=" * 104)
