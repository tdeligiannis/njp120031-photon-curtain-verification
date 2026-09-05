"""NJP-120031 v8 forward simulation -- core machinery.
Changes from v7: (i) time-parameterized kinematics -- make_psi/propagate/recon accept a
z->t map for the vertical free-fall geometry (transverse dynamics remain exactly free,
manuscript Sec. 2); constant-velocity mode retained for the Sec.-5.4 illustrative grid.
(ii) per-slit field components exposed (make_psi_parts) for the visibility mixture (m12).
(iii) HONEST Stage-1 gate (m13): the Sec.-5.4 grid statistics exclude x_c = 0, where the
true velocity vanishes by symmetry and relative error is undefined; x_c = 0 is reported
separately as an absolute-scale check. A second gate runs the identity in the actual
Configuration-E gravity geometry."""
import numpy as np, pickle
from scipy.interpolate import CubicSpline

def cum_integrate(integ, xc):
    """Cumulative integral of samples integ(xc) via cubic-spline antiderivative (O(h^4)).
    Replaces the v7 mixed Simpson/trapezoid rule, whose O(h^2) odd-step error at the
    fringe scale left a plane-independent ~2 um/s residual in the cumulative estimator
    (h/lambda_f is nearly constant across planes under the v8 kinematics)."""
    return CubicSpline(xc, integ).antiderivative()(xc)

hbar = 1.054571817e-34
M    = 1.4431606e-25          # 87Rb (Steck), standardized v8

def _tmap_or_const(vz, tmap):
    if tmap is not None: return tmap
    if vz is None: raise ValueError("need vz or tmap")
    return lambda z: np.asarray(z, dtype=float)/vz

def make_psi_parts(d, sig0, vz=None, tmap=None):
    """Per-slit analytic fields psi_a (slit at +d/2), psi_b (slit at -d/2).
    tmap: z -> elapsed time since the slit plane (gravity geometry); default z/vz."""
    tm = _tmap_or_const(vz, tmap)
    def s_t(z):
        t = tm(z)
        return sig0*(1.0 + 1j*hbar*t/(2*M*sig0**2))
    def psi_a(x, z):
        s = s_t(z); return np.exp(-(x - d/2)**2/(4*sig0*s))/np.sqrt(np.abs(s))
    def psi_b(x, z):
        s = s_t(z); return np.exp(-(x + d/2)**2/(4*sig0*s))/np.sqrt(np.abs(s))
    return psi_a, psi_b

def make_psi(d, sig0, vz=None, tmap=None):
    """Two-slit field, its x-derivative, and the true conditional velocity field."""
    tm = _tmap_or_const(vz, tmap)
    def s_t(z):
        t = tm(z)
        return sig0*(1.0 + 1j*hbar*t/(2*M*sig0**2))
    def psi(x, z):
        s = s_t(z)
        a = np.exp(-(x - d/2)**2/(4*sig0*s)); b = np.exp(-(x + d/2)**2/(4*sig0*s))
        return (a + b)/np.sqrt(np.abs(s))
    def dpsi(x, z):
        s = s_t(z)
        a = np.exp(-(x - d/2)**2/(4*sig0*s)); b = np.exp(-(x + d/2)**2/(4*sig0*s))
        return (-(x - d/2)/(2*sig0*s)*a - (x + d/2)/(2*sig0*s)*b)/np.sqrt(np.abs(s))
    def v_true(x, z):
        return hbar/M*np.imag(dpsi(x, z)/psi(x, z))
    return psi, dpsi, v_true

def propagate(field, x, dt):
    """Angular-spectrum free transverse propagation over an elapsed TIME dt.
    (v7 signature took (dz, vz); v8 callers pass dt = t(z2) - t(z1) directly.)"""
    N = x.size; dx = x[1] - x[0]
    k = 2*np.pi*np.fft.fftfreq(N, dx)
    return np.fft.ifft(np.fft.fft(field)*np.exp(-1j*hbar*k**2*dt/(2*M)))

def recon_noiseless(psi, x, z_j, x_c, Wk, L, vz=None, tmap=None, dphi=None):
    """v_x(x_c, z_j) from Eq.(13), plug-in denominator. Wk = transverse kernel width
    (v8 convention: W = w_eff). Time of flight t = t(L) - t(z_j) from the map."""
    tm = _tmap_or_const(vz, tmap)
    t = tm(L) - tm(z_j)
    ps = psi(x, z_j)
    if dphi is not None: ps = ps*np.exp(1j*dphi)
    f = np.exp(-(x - x_c)**2/Wk**2)
    G   = propagate(f*ps, x, t)
    psL = propagate(ps,   x, t)
    dx = x[1] - x[0]
    rho_w = np.sum(f*np.abs(ps)**2)*dx
    num = np.sum((x - x_c)*np.real(np.conj(psL)*G))*dx
    return num/(t*rho_w), rho_w, num

# =============================== Stage-1 gates ===============================
if __name__ == "__main__":
    # ---- Gate 1a: Sec.-5.4 illustrative grid (constant velocity, as in the text) ----
    # Geometry of manuscript Sec. 5.4: d = 5 um, sig0 = 0.6 um, W_kernel = 2 um,
    # v_z = 1 m/s, L = 0.30 m. HONEST statistics (m13): x_c = 0 excluded from
    # relative-error stats (v_true(0) = 0 by symmetry); reported as absolute check.
    d, sig0, Wk, vz, L = 5e-6, 0.6e-6, 2e-6, 1.0, 0.30
    psi, dpsi, v_true = make_psi(d, sig0, vz=vz)
    Xh, N = 2.0e-3, 2**16
    x = np.linspace(-Xh, Xh, N, endpoint=False)
    zs = [0.05, 0.10, 0.15, 0.20]; xcs = np.arange(-3, 4)*1e-6
    rows = []; errs = []; abs0 = []
    for z in zs:
        for xc in xcs:
            vr, _, _ = recon_noiseless(psi, x, z, xc, Wk, L, vz=vz)
            vd = v_true(np.array([xc]), z)[0]
            if abs(xc) < 1e-12:
                abs0.append((z, vd, vr)); rows.append((z, 0.0, vd, vr, np.nan)); continue
            e = abs(vr - vd)/abs(vd); errs.append(e)
            rows.append((z, xc*1e6, vd, vr, e*100))
    errs = np.array(errs)
    print("GATE 1a -- Sec.-5.4 grid, constant-v (24 points, x_c = 0 excluded honestly):")
    print(f"  mean rel err {errs.mean()*100:.2f}%   max {errs.max()*100:.2f}%   over {errs.size} points")
    a0 = max(max(abs(vd), abs(vr)) for _, vd, vr in abs0)
    vsc = np.median([abs(r[2]) for r in rows if np.isfinite(r[4])])
    print(f"  x_c = 0 absolute check: max |v| = {a0:.2e} m/s ({a0/vsc*100:.3f}% of median field scale)")
    g1a = errs.max() < 0.04
    # ---- Gate 1b: Configuration-E gravity geometry (t-parameterized map) ----
    import v8_config as C
    tmap = C.t_of_z
    psiE, _, v_trueE = make_psi(C.d, C.sig0, tmap=tmap)
    xE = np.linspace(-3e-3, 3e-3, 2**16, endpoint=False)
    cases = [(0.05, [-20e-6, -8e-6, 8e-6, 20e-6]),
             (0.30, [-80e-6, -30e-6, 30e-6, 80e-6]),
             (0.57, [-150e-6, -60e-6, 60e-6, 150e-6])]
    print("\nGATE 1b -- Configuration-E gravity map, Eq.(13) vs direct v_true:")
    eE = []
    for z, xcl in cases:
        for xc in xcl:
            vr, _, _ = recon_noiseless(psiE, xE, z, xc, C.W, C.L, tmap=tmap)
            vd = v_trueE(np.array([xc]), z)[0]
            e = abs(vr - vd)/abs(vd); eE.append(e)
        print(f"  z = {z*100:5.1f} cm : rel errs " +
              "  ".join(f"{abs(recon_noiseless(psiE,xE,z,xc,C.W,C.L,tmap=tmap)[0]-v_trueE(np.array([xc]),z)[0])/abs(v_trueE(np.array([xc]),z)[0])*100:.2f}%" for xc in xcl))
    eE = np.array(eE)
    print(f"  gravity-gate: mean {eE.mean()*100:.2f}%   max {eE.max()*100:.2f}%")
    g1b = eE.max() < 0.04
    print(f"\nSTAGE-1 GATES: 1a {'PASS' if g1a else 'FAIL'}   1b {'PASS' if g1b else 'FAIL'}")
    pickle.dump(dict(grid_rows=rows, grid_mean=errs.mean(), grid_max=errs.max(),
                     xc0_absmax=a0, gravity_mean=eE.mean(), gravity_max=eE.max(),
                     pass_1a=bool(g1a), pass_1b=bool(g1b)),
                open('sim/stage1_verification.pkl', 'wb'))
