"""NJP-120031 v8 simulation configuration -- single source of truth.
Everything here traces to the manuscript's verification scripts:
  kinematics  : verify_kinematics.py   (25/25, Gate G1)
  couplings   : verify_p2_couplings.py (68/68, Gate G2)  -- anchored at z_ref = 5 cm
  systematics : verify_p3_systematics.py (46/46, Gate G3) -- f-impurity, pairing
Convention (D3): w_spot = 2 um is the 1/e^2 intensity radius; the transverse
Gaussian KERNEL carries W = w_eff = w_spot/sqrt(2) (1/e intensity radius).
Per-plane scaling: g_a, N_gamma, Lambda_eff, P_sc  ~ tau_a(z) = w_spot/v_z(z);
phi_1 and U_0 are plane-independent; record noise SIG ~ tau_a^{-1/2}."""
import numpy as np

# ---------- fundamental / apparatus constants ----------
hbar = 1.054571817e-34
h    = 6.62607015e-34
M    = 1.4431606e-25          # 87Rb (Steck), standardized at Gate G0 (FLAG-6)
g    = 9.80665
v0   = 1.5                    # slit-plane longitudinal speed (primary spec, D1)
L    = 0.60                   # slit-to-screen drop
d    = 2.5e-6                 # slit separation
sig0 = 0.6e-6                 # slit Gaussian sigma
w_spot = 2.0e-6               # 1/e^2 intensity radius (laser spec)
W      = w_spot/np.sqrt(2.0)  # 1/e intensity radius: ALL transverse shape factors
lamD2  = 780.241209686e-9
vrec   = h/(M*lamD2)          # 5.884 mm/s

# ---------- kinematics (Gate G1) ----------
def vz_of(z):  return np.sqrt(v0**2 + 2.0*g*np.asarray(z, dtype=float))
def t_of_z(z): return (vz_of(z) - v0)/g
T       = float(t_of_z(L))            # 0.228831 s
T_FALL  = v0/g                        # 0.152957 s release->slit (v0 from rest)
FRINGE  = h*T/(M*d)                   # 420.3 um accelerated-fringe period
def lam_f(z):  return h*t_of_z(z)/(M*d)   # local fringe wavelength at plane z
def tau_of(z): return w_spot/vz_of(z)
Z_REF   = 0.05
TAU_REF = float(tau_of(Z_REF))        # 1.1127 us

# ---------- couplings anchored at z_ref (Gate G2, Table 2) ----------
GA_M_REF, GA_C_REF = 3673.0, -3679.0
NG_M_REF, NG_C_REF = 1.40e9, 3.76e9
LEFF_REF, PSC_REF  = 8.65e-3, 0.0144
def _s(z): return tau_of(z)/TAU_REF
def ga_m_of(z): return GA_M_REF*_s(z)
def ga_c_of(z): return GA_C_REF*_s(z)
def Ng_m_of(z): return NG_M_REF*_s(z)
def Ng_c_of(z): return NG_C_REF*_s(z)
def Leff_of(z): return LEFF_REF*_s(z)
def V_of(z):    return np.exp(-Leff_of(z))       # fringe visibility (m12)
def Psc_of(z):  return PSC_REF*_s(z)
def _sig_arm(Ng, ga): return np.sqrt(Ng)/(np.sqrt(2.0)*abs(ga))
def SIG_of(z):
    """Combined per-shot dimensionless record noise (units hbar*U0/W); ~ tau_a^{-1/2}."""
    sm = _sig_arm(Ng_m_of(z), ga_m_of(z)); sc = _sig_arm(Ng_c_of(z), ga_c_of(z))
    return 1.0/np.sqrt(sm**-2 + sc**-2)

# ---------- f-impurity term (Gate G3, sec:magnetic; forward-model obligation) ----------
F_IMP  = 0.05                 # Zeeman-impurity fraction at spec
A_MAG  = 3.213e-5             # m/s^2 at spec grad 0.1 mG/cm, |gF mF| = 1/2
# conservative screen displacement: post-slit (0.5 a T^2 = 0.84 um) + pre-slit
# velocity accrual carried through (a t_fall T = 1.12 um); split +-mF, f/2 each
DSCR_IMP = 0.5*A_MAG*T**2 + A_MAG*T_FALL*T       # 1.97 um total

# ---------- record-atom pairing ceiling (Gate G3, sec:detection; CONFLICT-2) ----------
DU_PAIR  = 5.0e-3             # longitudinal Raman selection width (m/s)
DU_THERM = 1.5e-2             # unselected longitudinal thermal width
FLUX_LO, FLUX_HI = 3.2e6, 8.6e7   # source window, atoms/day (pre-selection)
def pair_coef(z): return (vz_of(L) - vz_of(z))/(g*vz_of(z))   # s per (m/s)
def ceiling_day(z, du=DU_PAIR): return 86400.0/(2.0*np.e*pair_coef(z)*du)
def eff_window(z):
    """Per-plane usable window [lo, hi] atoms/day after pairing veto, optimizing
    the selection width du <= DU_THERM independently at each flux edge."""
    out = []
    for edge in (FLUX_LO, FLUX_HI):
        du_star = np.sqrt(86400.0/(2.0*np.e*pair_coef(z))*DU_THERM/edge)
        du_star = min(du_star, DU_THERM)
        out.append(min(edge*du_star/DU_THERM, ceiling_day(z, du_star)))
    return tuple(out)

# ---------- self-checks against the manuscript verification scripts ----------
def _check():
    a = lambda x, r, tol=5e-3: abs(x - r) <= tol*abs(r)
    assert a(T, 0.228831), T
    assert a(FRINGE, 420.3e-6), FRINGE
    assert a(vz_of(L), 3.7441), vz_of(L)
    assert a(TAU_REF, 1.1127e-6), TAU_REF
    assert a(GA_M_REF/NG_M_REF, 2.623e-6), "phi1_m"      # plane-independent
    assert a(abs(GA_C_REF)/NG_C_REF, 9.785e-7), "phi1_c"
    assert a(SIG_of(Z_REF), 6.147), SIG_of(Z_REF)
    assert a(SIG_of(0.57), 6.147*np.sqrt(TAU_REF/tau_of(0.57))), "SIG scaling"
    assert a(vrec, 5.884e-3), vrec
    assert a(pair_coef(Z_REF), 0.11044), pair_coef(Z_REF)
    assert a(ceiling_day(Z_REF), 2.88e7, 1e-2), ceiling_day(Z_REF)
    assert a(V_of(Z_REF), 0.99139, 1e-4), V_of(Z_REF)
    assert a(DSCR_IMP, 1.97e-6, 2e-2), DSCR_IMP
_check()

if __name__ == "__main__":
    print("v8_config self-checks PASS")
    print(f"  T = {T:.6f} s   FRINGE = {FRINGE*1e6:.1f} um   v_z(L) = {vz_of(L):.4f} m/s")
    zj = np.linspace(0.05, 0.57, 20)
    print(f"  SIG(z): {SIG_of(0.05):.3f} (z_ref) -> {SIG_of(0.57):.3f} (last)")
    print(f"  V(z):   {V_of(0.05):.4f} -> {V_of(0.57):.4f}    Psc(z): {Psc_of(0.05):.4f} -> {Psc_of(0.57):.4f}")
    print(f"  pairing eff_window z=5cm: {eff_window(0.05)[0]:.2e} - {eff_window(0.05)[1]:.2e} /day"
          f"   (unconstrained window {FLUX_LO:.1e} - {FLUX_HI:.1e})")
    print(f"  eff_window z=30cm: {eff_window(0.2963)[0]:.2e} - {eff_window(0.2963)[1]:.2e} /day")
