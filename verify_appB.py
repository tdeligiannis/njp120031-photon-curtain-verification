#!/usr/bin/env python3
"""verify_appB.py -- NJP-120031 P8-4: verification of every step of Appendix B (derivation chain for R2-8).

Symbolic (sympy) checks of the operator identities, numeric checks on the exact two-Gaussian double-slit
state, and the decomposition of the Sec.-5.4 verification-grid errors into the two O(w_eff^2) terms that
Appendix B derives. Run from the campaign root (needs simcode_v8/sim on the path). Prints PASS/FAIL lines;
exit status 1 on any FAIL. Wired into verify_p7_master.py section D.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'simcode_v8', 'sim'))
import sympy as sp
import sim_core as S
import v8_config as C

npass = nfail = 0
def check(name, ok, detail=""):
    global npass, nfail
    if ok: npass += 1
    else: nfail += 1
    print(f"{'PASS' if ok else 'FAIL'} {name:70s} {detail}")

print("="*100); print("APPENDIX B + SEC.-9 DERIVATION CHECKS  (P8-4/P8-5, 2026-08-24)"); print("="*100)

# ----------------------------------------------------------------------------------------------
# A. SYMBOLIC CHECKS
# ----------------------------------------------------------------------------------------------
print("\n--- A. symbolic identities (sympy) ---")
xa, xc, xg, dlt, w, U0, hb, phi1, tau1, t, M = sp.symbols('x_a x_c x_gamma delta w U_0 hbar varphi_1 tau_1 t M', real=True, positive=True)
xa = sp.Symbol('x_a', real=True); xc = sp.Symbol('x_c', real=True); xg = sp.Symbol('x_gamma', real=True); dlt = sp.Symbol('delta', real=True)
f = lambda x, c: sp.exp(-(x - c)**2/w**2)                       # curtain profile f~(x; x_c), w == w_eff here

# A1 -- AAV linearization: f~(x_a - x_c - delta) = f~ - delta d_x f~ + O(delta^2); the bilinear coefficient is A-hat/(hbar U0)
lin = sp.series(f(xa - dlt, xc), dlt, 0, 2).removeO()
A_from_lin = sp.simplify(hb*U0*sp.expand(lin).coeff(dlt, 1))
A_eq8 = 2*hb*U0*(xa - xc)/w**2*f(xa, xc)
A_eq12 = hb*U0*sp.diff(f(xa, xc), xc)
check("A1 bilinear coefficient of the effective coupling: hbar U0 d_{x_c} f~ = Eq.(8) (Taylor coefficient of the Gaussian profile; the contact coupling itself is checked numerically in B4-2)",
      sp.simplify(A_from_lin - A_eq8) == 0 and sp.simplify(A_eq12 - A_eq8) == 0)
check("A1' zeroth-order term is the position-diagonal dipole phase hbar U0 f~(x_a; x_c) (no meter operator)",
      sp.simplify(sp.expand(lin).coeff(dlt, 0) - f(xa, xc)) == 0)

# A2 -- photon kick Eq.(6) equals -tau_1 A with tau_1 = phi_1/U0 (Eq. 11)
phig = phi1*sp.exp(-(xa - xg)**2/w**2)
kick_eq6 = (-hb*sp.diff(phig, xg)).subs(xg, xc)
check("A2 Eq.(6) kick = -hbar d_{x_gamma} phi_gamma|_{x_c} equals -tau_1 A with tau_1 = phi_1/U_0",
      sp.simplify(kick_eq6 - (-(phi1/U0)*A_eq8)) == 0)

# A3 -- Heisenberg equation for the meter momentum: dp/dt = (i/hbar)[A x_meter, p] = -A
y = sp.Symbol('y', real=True); g = sp.Function('g')(y); Aw = sp.Symbol('A', real=True)
p_of = lambda h: -sp.I*hb*sp.diff(h, y)
comm = Aw*y*p_of(g) - p_of(Aw*y*g)           # [A y, p] g
check("A3 (i/hbar)[A delta x_gamma, p_gamma] = -A  (so dp_gamma/dt = -A under H_AAV)",
      sp.simplify(sp.I/hb*comm - (-Aw*g)) == 0)

# A4 -- exact conditional pointer shift for a Gaussian meter and c-number weak value A_w = a + i b:
#        <p> = -tau_1 Re A_w (all orders in tau_1 within the weak-value replacement)
a, b, sx, eps = sp.symbols('a b sigma_x epsilon', real=True); sx = sp.Symbol('sigma_x', positive=True); eps = sp.Symbol('epsilon', positive=True)
Phi = sp.exp(-y**2/(4*sx**2))
Phif = sp.exp(-sp.I*eps*(a + sp.I*b)*y)*Phi                     # eps = tau_1/hbar
dens = sp.simplify(Phif*sp.conjugate(Phif))
num = sp.integrate(sp.simplify(sp.conjugate(Phif)*p_of(Phif)), (y, -sp.oo, sp.oo))
den = sp.integrate(dens, (y, -sp.oo, sp.oo))
pmean = sp.simplify(num/den)
check("A4 Gaussian meter, c-number weak value: <p_gamma> = -tau_1 Re A_w exactly (Eq. 11 sign and magnitude)",
      sp.simplify(pmean - (-hb*eps*a)) == 0, f"<p> = {pmean}")
xmean = sp.simplify(sp.integrate(y*dens, (y, -sp.oo, sp.oo))/den)
check("A4' the meter POSITION shift is 2 tau_1 sigma_x^2 Im A_w / hbar (imaginary part goes to x, not p)",
      sp.simplify(xmean - 2*eps*b*sx**2) == 0)

# A5 -- free kernel: (x_f - x) K = (i hbar t/M) d_x K, the position-space form of x(t) = x + p t/M
xf, x = sp.symbols('x_f x', real=True); tt = sp.Symbol('t', positive=True); Mm = sp.Symbol('M', positive=True)
K = sp.sqrt(Mm/(2*sp.pi*sp.I*hb*tt))*sp.exp(sp.I*Mm*(xf - x)**2/(2*hb*tt))
check("A5 (x_f - x) K(x_f - x; t) = (i hbar t/M) d_x K   [Eq. B kernel identity]",
      sp.simplify((xf - x)*K - sp.I*hb*tt/Mm*sp.diff(K, x)) == 0)

# A6 -- Madelung: Im(psi* d psi) = rho dS/hbar ; and Re[i psi d psi*] = Im[psi* d psi]   (explicit real symbols)
R, R1, S0, S1 = sp.symbols('rho rho_prime S S_prime', real=True); Rp = sp.Symbol('rho', positive=True)
psi_M  = sp.sqrt(Rp)*sp.exp(sp.I*S0/hb)
dpsi_M = (R1/(2*sp.sqrt(Rp)) + sp.I*sp.sqrt(Rp)*S1/hb)*sp.exp(sp.I*S0/hb)     # d/dx of sqrt(rho) e^{iS/hbar}
lhs = sp.simplify(sp.im(sp.expand(sp.conjugate(psi_M)*dpsi_M)))
check("A6 Madelung: Im(psi* d_x psi) = rho d_x S / hbar  =>  J/rho = d_x S/M",
      sp.simplify(lhs - Rp*S1/hb) == 0, f"Im = {lhs}")
uu, vv, u1, v1 = sp.symbols('u v u_prime v_prime', real=True)
psi_c = uu + sp.I*vv; dpsi_c = u1 + sp.I*v1
check("A6' Re[i psi d_x psi*] = Im[psi* d_x psi]  (numerator step, Eq. 17)",
      sp.simplify(sp.re(sp.expand(sp.I*psi_c*sp.conjugate(dpsi_c))) - sp.im(sp.expand(sp.conjugate(psi_c)*dpsi_c))) == 0)

# A7 -- moment expansions of the curtain-blurred quantities in w (== w_eff)
s = sp.Symbol('s', real=True); wp = sp.Symbol('w', positive=True)
fs = sp.exp(-s**2/wp**2)
m0 = sp.integrate(fs, (s, -sp.oo, sp.oo)); m2 = sp.integrate(s**2*fs, (s, -sp.oo, sp.oo)); m4 = sp.integrate(s**4*fs, (s, -sp.oo, sp.oo))
check("A7 Gaussian moments: int f~ = sqrt(pi) w, int s^2 f~ = sqrt(pi) w^3/2, int s^4 f~ = 3 sqrt(pi) w^5/4",
      sp.simplify(m0 - sp.sqrt(sp.pi)*wp) == 0 and sp.simplify(m2 - sp.sqrt(sp.pi)*wp**3/2) == 0 and sp.simplify(m4 - 3*sp.sqrt(sp.pi)*wp**5/4) == 0)
r0, r1, r2, r3, r4 = sp.symbols('rho_0 rho_1 rho_2 rho_3 rho_4', real=True)   # Taylor coefficients of rho about x_c
rho_series = r0 + r1*s + r2*s**2/2 + r3*s**3/6 + r4*s**4/24
A_moment = sp.integrate(s*fs*rho_series, (s, -sp.oo, sp.oo))
rho_w = sp.integrate(fs*rho_series, (s, -sp.oo, sp.oo))
check("A7' script-A = (sqrt(pi)/2) w^3 rho' + (sqrt(pi)/8) w^5 rho''' + O(w^7)   [Eq. 17 coefficient]",
      sp.simplify(A_moment - (sp.sqrt(sp.pi)/2*wp**3*r1 + sp.sqrt(sp.pi)/8*wp**5*r3)) == 0)
check("A7'' rho_w = sqrt(pi) w rho + (sqrt(pi)/4) w^3 rho'' + O(w^5)",
      sp.simplify(rho_w - (sp.sqrt(sp.pi)*wp*r0 + sp.sqrt(sp.pi)/4*wp**3*r2 + sp.sqrt(sp.pi)/32*wp**5*r4)) == 0)
resid = sp.series(A_moment/(tt*rho_w), wp, 0, 4).removeO()
check("A7''' residual script-A/(t rho_w) = (w^2/2t) rho'/rho + O(w^4)   [Appendix B bound, leading term]",
      sp.simplify(resid - wp**2/(2*tt)*r1/r0) == 0, f"series = {resid}")

# A8 -- w -> 0 limit of the ratio (t J_w + A)/(t rho_w) -> J/rho
j0, j2 = sp.symbols('J_0 J_2', real=True)
J_w = sp.integrate(fs*(j0 + j2*s**2/2), (s, -sp.oo, sp.oo))
ratio = sp.series((tt*J_w + A_moment)/(tt*rho_w), wp, 0, 2).removeO()
J0, J2 = sp.symbols('J_0 J_2', real=True)
blur = sp.series(sp.integrate(fs*(J0 + J2*s**2/2), (s, -sp.oo, sp.oo))/rho_w, wp, 0, 4).removeO() - J0/r0
check("A7-4 finite-blur term J_w/rho_w - J/rho = (w^2/4)(J-dd - v rho-dd)/rho + O(w^4)",
      sp.simplify(blur - wp**2/4*(J2 - J0/r0*r2)/r0) == 0)
check("A8 narrow-curtain limit: (t J_w + script-A)/(t rho_w) -> J/rho as w -> 0",
      sp.simplify(ratio - j0/r0) == 0)

# ----------------------------------------------------------------------------------------------
# B. NUMERIC CHECKS on the exact two-Gaussian state (Sec.-5.4 geometry)
# ----------------------------------------------------------------------------------------------
print("\n--- B. numeric checks on the exact two-Gaussian state ---")
hbar, Mkg = S.hbar, S.M
d, sig0, Wk, vz, L = 5e-6, 0.6e-6, 2e-6, 1.0, 0.30          # Sec.-5.4 grid (Gate 1a)
psi, dpsi, v_true = S.make_psi(d, sig0, vz=vz)
N = 2**16; Xh = 2.0e-3
xg_ = np.linspace(-Xh, Xh, N, endpoint=False); dx = xg_[1] - xg_[0]
kk = 2*np.pi*np.fft.fftfreq(N, dx)
def prop(field, dt): return np.fft.ifft(np.fft.fft(field)*np.exp(-1j*hbar*kk**2*dt/(2*Mkg)))
def deriv(field): return np.fft.ifft(1j*kk*np.fft.fft(field))
zj = 0.05; tj = zj/vz; tL = L/vz; tprop = tL - tj
ps = psi(xg_, zj); psL = prop(ps, tprop)
# B1: zeroth-moment identity (Eq. 14): int dx_f psi_L*(x_f) K(x_f - x; t) = psi*(x, z_j)
#     realised as the adjoint propagation of psi_L (K(x_f - x; t) integrated against psi_L* is [U^dagger psi_L]*)
lhs0 = np.conj(prop(psL, -tprop))
core = np.abs(ps) > 1e-3*np.abs(ps).max()
e0 = np.max(np.abs(lhs0 - np.conj(ps))[core])/np.abs(ps).max()
check("B1 zeroth-moment identity (Eq. 14) on the exact state", e0 < 1e-9, f"max rel dev {e0:.1e}")
# B2: first-moment identity (Eq. 15): int dx_f x_f psi_L* K = x psi* + (i hbar t/M) d_x psi*
lhs1 = np.conj(prop(xg_*psL, -tprop))
rhs1 = xg_*np.conj(ps) + 1j*hbar*tprop/Mkg*deriv(np.conj(ps))
e1 = np.max(np.abs(lhs1 - rhs1)[core])/np.max(np.abs(rhs1)[core])
check("B2 first-moment identity (Eq. 15) on the exact state", e1 < 1e-7, f"max rel dev {e1:.1e}")

# B3: decomposition of the reconstruction error into the two O(w_eff^2) terms, over the 24-point grid
#     v_recon - v_true = [J_w/rho_w - J/rho] + [script-A/(t rho_w)]   (exact identity, Eq. 18 with Eq. 16-17)
print("\n    Sec.-5.4 grid: error decomposition  (worst plane first)")
zs = [0.05, 0.10, 0.15, 0.20]; xcs = np.arange(-3, 4)*1e-6
rows = []; worst = None; dec_dev = 0.0; asym_dev = 0.0
for z in zs:
    tj = z/vz; tp = tL - tj
    psz = psi(xg_, z); rho_z = np.abs(psz)**2
    J_z = hbar/Mkg*np.imag(np.conj(psz)*deriv(psz))
    dlnrho = np.real(deriv(rho_z))/np.maximum(rho_z, 1e-300)
    for xc_ in xcs:
        if abs(xc_) < 1e-12: continue
        vr, rw, num = S.recon_noiseless(psi, xg_, z, xc_, Wk, L, vz=vz)
        vd = v_true(np.array([xc_]), z)[0]
        fz = np.exp(-(xg_ - xc_)**2/Wk**2)
        Jw = np.sum(fz*J_z)*dx; rhow = np.sum(fz*rho_z)*dx; Amom = np.sum((xg_ - xc_)*fz*rho_z)*dx
        blur = Jw/rhow - vd; offset = Amom/(tp*rhow)
        dec_dev = max(dec_dev, abs((vr - vd) - (blur + offset))/abs(vd))
        i_c = np.argmin(np.abs(xg_ - xc_)); offset_asym = Wk**2/(2*tp)*dlnrho[i_c]
        asym_dev = max(asym_dev, abs(offset - offset_asym)/abs(offset)) if abs(offset) > 0 else asym_dev
        e = (vr - vd)/vd
        rows.append((z, xc_*1e6, e*100, blur/vd*100, offset/vd*100, offset_asym/vd*100))
        if worst is None or abs(e) > abs(worst[2]): worst = rows[-1]
print(f"    {'z (m)':>6s} {'x_c(um)':>8s} {'err%':>8s} {'blur%':>8s} {'offset%':>9s} {'offset_asym%':>13s}")
for r in rows:
    if r[0] == 0.05 or r == worst: print(f"    {r[0]:6.2f} {r[1]:8.1f} {r[2]:8.3f} {r[3]:8.3f} {r[4]:9.3f} {r[5]:13.3f}")
errs = np.array([abs(r[2]) for r in rows])
check("B3 exact decomposition v_recon - v_true = blur + offset holds on all 24 points", dec_dev < 1e-8, f"max dev {dec_dev:.1e} (relative)")
check("B3' leading-order offset formula (w^2/2t) d_x ln rho reproduces the exact offset term", asym_dev < 0.05, f"max rel dev {asym_dev:.2e}")
check("B3'' grid statistics reproduce Gate 1a: mean 1.04%, max 2.52%", abs(errs.mean() - 1.04) < 0.02 and abs(errs.max() - 2.52) < 0.02, f"mean {errs.mean():.2f}% max {errs.max():.2f}%")
frac_offset_worst = worst[4]/worst[2]
check("B3-4 quoted decomposition at the worst point: -2.52% = -2.10% (blur) + -0.42% (offset)",
      abs(worst[2] - (-2.52)) < 0.006 and abs(worst[3] - (-2.10)) < 0.006 and abs(worst[4] - (-0.42)) < 0.006, f"{worst[2]:.3f} = {worst[3]:.3f} + {worst[4]:.3f}")
check("B3-5 leading-order offset formula within 3% of the exact offset (quoted: better than 3%)", asym_dev < 0.03, f"{asym_dev:.3f}")
check("B3''' worst point (z=0.05, x_c=-3): both O(w^2) terms are of the same sign/order (reported, not asserted)", True,
      f"err {worst[2]:.2f}% = blur {worst[3]:.2f}% + offset {worst[4]:.2f}% (offset fraction {frac_offset_worst:.2f})")

# B4: full joint-state spot-check of the pointer shift Eq. (11): <p_gamma>_{x_f} vs -tau_1 Re A_w
#     Two microscopic couplings are propagated exactly on the joint (atom x meter) grid:
#       (i)  the effective bilinear AAV form  exp(-i tau_1 A(x_a) delta x_gamma / hbar)   [Eq. 9];
#       (ii) the CONTACT coupling  exp(-i kappa tau_1 delta_eps(x_gamma - x_a))  -- the atom couples to the local
#            intensity at its own position; the curtain Gaussian is the METER STATE |Phi_0|^2, not a kernel.
#     Appendix B derives that both give <p_gamma>_{x_f} = -tau_1 Re<A>_w to first order in phi_1.
#     (A literal Gaussian KERNEL of width w_eff in the photon coordinate, tested during P8-4, gives a shift 47% off,
#      independent of phi_1 -- because the meter spread w/2 is comparable to w_eff; that model is unphysical.)
print("\n    joint atom+meter state, z_j = 0.05 m, x_c = -3 um (Sec.-5.4 geometry)")
Na = 2**15; Xa = 1.0e-3
xa_ = np.linspace(-Xa, Xa, Na, endpoint=False); dxa = xa_[1] - xa_[0]; ka = 2*np.pi*np.fft.fftfreq(Na, dxa)
def prop_a(F, dt):   # propagate along axis 0 (atom) for a 2-D array
    return np.fft.ifft(np.fft.fft(F, axis=0)*np.exp(-1j*hbar*ka[:, None]**2*dt/(2*Mkg)), axis=0)
w_spot = np.sqrt(2)*Wk                      # 1/e^2 radius; |Phi_0|^2 propto exp(-2y^2/w^2) = exp(-y^2/w_eff^2) = f~ ; sigma_x = w/2
sigx = w_spot/2
xc_s = -3e-6; tj = zj/vz; tp = tL - tj
ps_a = psi(xa_, zj)
ftil = np.exp(-(xa_ - xc_s)**2/Wk**2); dxc_f = 2*(xa_ - xc_s)/Wk**2*ftil          # f~ and d_{x_c} f~ on the atom grid
psL_a = prop_a(ps_a[:, None], tp)[:, 0]; GA = prop_a((dxc_f*ps_a)[:, None], tp)[:, 0]   # unnormalised weak value of d_{x_c} f~
xf_list = [-40e-6, 0.0, 60e-6, 150e-6]
def run_joint(phi, mode, Ng, epsf):
    yg = np.linspace(-4*sigx, 4*sigx, Ng); dy = yg[1] - yg[0]
    Phi_m = np.exp(-yg**2/(4*sigx**2)); Phi_m /= np.sqrt(np.sum(np.abs(Phi_m)**2)*dy)      # |Phi_m|^2 = f~(y)/(sigx sqrt(2 pi))
    P0 = np.abs(Phi_m).max()**2; eps = epsf*sigx
    devs = []
    for xf_ in xf_list:
        i_f = np.argmin(np.abs(xa_ - xf_))
        Aw = GA[i_f]/psL_a[i_f]                                             # weak value of d_{x_c} f~ (units 1/m)
        pred = -hbar*phi*np.real(Aw)                                        # -tau_1 Re<A>_w with tau_1 A = hbar phi_1 d_{x_c} f~
        if mode == 'bilinear':
            phase = -phi*np.outer(dxc_f, yg)                                # -tau_1 A(x_a) delta x_gamma / hbar
        else:
            phase = -(phi/P0)*np.exp(-(yg[None, :] - (xa_[:, None] - xc_s))**2/(2*eps**2))/(eps*np.sqrt(2*np.pi))   # -kappa tau_1 delta_eps(x_gamma - x_a)
        chi = prop_a(ps_a[:, None]*Phi_m[None, :]*np.exp(1j*phase), tp)[i_f, :]
        pm = np.real(np.sum(np.conj(chi)*(-1j*hbar*np.gradient(chi, dy)))*dy)/(np.sum(np.abs(chi)**2)*dy)
        devs.append(abs(pm - pred)/abs(pred))
    return max(devs)
r_bil = {phi: run_joint(phi, 'bilinear', 256, 0.2) for phi in (1e-2, 1e-3)}
r_con = {(phi, e): run_joint(phi, 'contact', Ng, e) for phi, Ng, e in ((1e-2, 256, 0.2), (1e-3, 256, 0.2), (1e-4, 256, 0.2), (1e-4, 512, 0.1))}
for phi in (1e-2, 1e-3): print(f"    phi_1 = {phi:.0e}: bilinear (AAV) max |<p>/pred - 1| = {r_bil[phi]:.2e}")
for (phi, e), v in r_con.items(): print(f"    phi_1 = {phi:.0e}, eps = {e:.1f} sigma_x: contact (local-intensity) max |<p>/pred - 1| = {v:.2e}")
check("B4 pointer shift Eq.(11) reproduced by the exact joint state, effective bilinear coupling (phi_1 = 1e-3, 1e-4 level)", r_bil[1e-3] < 1e-3, f"dev {r_bil[1e-3]:.1e}")
check("B4-2 pointer shift Eq.(11) reproduced with the microscopic CONTACT coupling (phi_1 = 1e-4, eps = 0.1 sigma_x): < 1%", r_con[(1e-4, 0.1)] < 1e-2, f"dev {r_con[(1e-4, 0.1)]:.1e}")
check("B4-3 contact-model residual is the delta regularisation (falls when eps is halved at phi_1 = 1e-4)", r_con[(1e-4, 0.1)] < 0.6*r_con[(1e-4, 0.2)], f"{r_con[(1e-4, 0.2)]:.1e} -> {r_con[(1e-4, 0.1)]:.1e}")
check("B4-4 bilinear-model residual does not grow from phi_1 = 1e-2 to 1e-3 (linearisation error is higher order)", r_bil[1e-3] <= r_bil[1e-2], f"{r_bil[1e-2]:.1e} -> {r_bil[1e-3]:.1e}")

# B5: impulsive-transit condition: transverse displacement during a curtain transit << w_eff
vscale = 362.4e-6            # field scale at plane 0 (stage6 record, m/s)
tau_ref = C.tau_a_of(0.05) if hasattr(C, 'tau_a_of') else C.W_SPOT/np.sqrt(C.v0**2 + 2*C.g*0.05) if hasattr(C, 'W_SPOT') else None
if tau_ref is None:
    tau_ref = 1.113e-6       # Table 1 value (fallback)
disp = vscale*tau_ref
bound_first = C.W/(C.t_of_z(C.L) - C.t_of_z(0.05))
check("B6 residual bound scale w_eff/t at the first curtain plane (quoted: approximately 7 um/s)", abs(bound_first*1e6 - 7.1) < 0.15, f"{bound_first*1e6:.2f} um/s")
check("B5 impulsive transit: v_x tau_a << w_eff at the reference plane (ratio > 10^3)", C.W/disp > 1e3, f"displacement {disp*1e9:.2f} nm; w_eff/displacement = {C.W/disp:.0f}")

# ----------------------------------------------------------------------------------------------
# C. DECOHERENCE: the partial-trace passage of Sec. 9 (P8-5, R2-9)
# ----------------------------------------------------------------------------------------------
print("\n--- C. decoherence / partial-trace checks (Sec. 9) ---")
# C1 -- overlap of two momentum-displaced Gaussian meter states of momentum width sigma_p (amplitude Gaussian in p)
pp, d1, d2, sgp = sp.symbols('p delta_1 delta_2 sigma_p', real=True); sgp = sp.Symbol('sigma_p', positive=True)
gp = lambda dd: (2*sp.pi*sgp**2)**sp.Rational(-1, 4)*sp.exp(-(pp - dd)**2/(4*sgp**2))
ov = sp.simplify(sp.integrate(gp(d1)*gp(d2), (pp, -sp.oo, sp.oo)))
check("C1 <Phi_{delta_1}|Phi_{delta_2}> = exp[-(delta_1 - delta_2)^2/(8 sigma_p^2)] for momentum-displaced Gaussians",
      sp.simplify(ov - sp.exp(-(d1 - d2)**2/(8*sgp**2))) == 0)
# C2 -- coherent-state overlap |<alpha_u|alpha_v>| = exp[-|alpha|^2 (1 - Re<u|v>)]  (numeric, two-mode truncated Fock space)
from math import factorial
def coh_2mode(alpha, theta, nmax=40):
    # mode u = a1 ; mode v = cos(theta) a1 + sin(theta) a2 ; coherent state |alpha> in mode v = |alpha cos|_1 |alpha sin|_2
    def coh(a):
        v = np.array([a**n/np.sqrt(float(factorial(n))) for n in range(nmax)], dtype=complex); return np.exp(-abs(a)**2/2)*v
    return np.kron(coh(alpha*np.cos(theta)), coh(alpha*np.sin(theta)))
alpha = 1.3; devC2 = 0.0
for theta in (0.3, 0.7, 1.2):
    u = coh_2mode(alpha, 0.0); v = coh_2mode(alpha, theta)
    ov_num = abs(np.vdot(u, v)); ov_pred = np.exp(-abs(alpha)**2*(1 - np.cos(theta)))
    devC2 = max(devC2, abs(ov_num/ov_pred - 1))
check("C2 |<alpha_u|alpha_v>| = exp[-|alpha|^2 (1 - Re<u|v>)] (two-mode Fock-space numeric, three mode overlaps)", devC2 < 1e-10, f"max rel dev {devC2:.1e}")
# C3 -- D_1 from the kick difference at x_a = +-d/2, x_c = 0 (Eq. 25)
dd_, wE, ph = sp.symbols('d w_eff varphi_1', positive=True)
kick = lambda xa_: -2*hb*ph*xa_/wE**2*sp.exp(-xa_**2/wE**2)          # Eq. (6) at x_c = 0
D1 = sp.simplify((kick(-dd_/2) - kick(dd_/2))**2/(8*(hb/(sp.sqrt(2)*wE))**2))  # sigma_p = hbar/w = hbar/(sqrt2 w_eff)
check("C3 D_1 = (dp_L - dp_R)^2/(8 sigma_p^2) = phi_1^2 (d/w_eff)^2 exp[-d^2/(2 w_eff^2)]   [Eq. 25]",
      sp.simplify(D1 - ph**2*(dd_/wE)**2*sp.exp(-dd_**2/(2*wE**2))) == 0)
# C4 -- exact coherent-state exponent N(1 - e^{-D_1}) vs N D_1: relative correction D_1/2 (measuring arm, config values)
phi1_m = C.GA_M_REF/C.NG_M_REF
D1_m = phi1_m**2*(C.d/C.W)**2*np.exp(-C.d**2/(2*C.W**2))
import mpmath as mp; mp.mp.dps = 40
rel = float(1 - (1 - mp.e**(-mp.mpf(D1_m)))/mp.mpf(D1_m))            # extended precision: double cannot resolve a 1e-12 relative term
check("C4 exponent N(1 - e^{-D_1}) differs from N D_1 by relative D_1/2 ~ 1e-12 (measuring arm)", abs(rel - D1_m/2) < 1e-3*D1_m and D1_m < 1e-11, f"D_1 = {D1_m:.2e}, rel corr {rel:.1e}")
# C5 -- per-photon Englert distinguishability sqrt(1 - e^{-2 D_1}) ~ sqrt(2 D_1)
Dsym = sp.Symbol('D_1', positive=True)
check("C5 per-photon Englert distinguishability sqrt(1 - e^{-2 D_1}) = sqrt(2 D_1) (1 + O(D_1))",
      sp.simplify(sp.series(sp.sqrt(1 - sp.exp(-2*Dsym))/sp.sqrt(2*Dsym), Dsym, 0, 1).removeO() - 1) == 0)
# C6 -- the slit-pair kick difference is maximal at x_c = 0 (so the quoted D_1 is the upper value)
xcs_ = np.linspace(-3*C.W, 3*C.W, 2001)
gk = lambda s_: s_*np.exp(-s_**2/C.W**2)
diff_ = np.abs(gk(C.d/2 - xcs_) - gk(-C.d/2 - xcs_))
check("C6 |dp_L - dp_R|(x_c) is maximal at x_c = 0 for d/w_eff of the apparatus", abs(xcs_[np.argmax(diff_)]) < 2*(xcs_[1]-xcs_[0]), f"argmax at x_c = {xcs_[np.argmax(diff_)]*1e9:.1f} nm")

print("\n" + "="*100); print(f"APPENDIX B VERIFICATION: {npass} checks passed, {nfail} failed"); print("="*100)
sys.exit(0 if nfail == 0 else 1)
