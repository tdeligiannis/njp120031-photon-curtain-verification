# NJP-120031 v8 simulation chain -- LOCAL RUN PACKAGE
Fabianiak & Deligiannis. Built 2026-07-04 from the v7 supplementary bundle for the v8
correction campaign (Gates G1-G3 closed; this run feeds P4).

## What changed from v7 (all edits asserted-unique; builder: build_simcode_v8.py)
1. KINEMATICS (G1): time-parameterized z<->t map everywhere (sim/v8_config.py:
   vz(z), t(z), T = 0.228831 s, fringe 420.3 um). Trajectory ODEs use dx/dz = v_x/v_z(z).
2. CONVENTION (G2/D3): transverse kernel width W = w_eff = w/sqrt(2) = 1.414 um in all
   shape factors, /W moment normalization, per-plane couplings anchored at z_ref = 5 cm
   (ga 3673/-3679, Ng 1.40e9/3.76e9, SIG(z_ref) = 6.15 scaling as tau_a^{-1/2}).
3. m12: screen sampling from the visibility mixture rho_V = V(z) rho_coh + (1-V) rho_inc,
   V(z) = exp(-Lambda_eff(z)); per-slit fields via make_psi_parts.
4. SE background: per-plane P_sc(z) times the recoil-smeared coherent density
   (box +-v_rec (T - t(z)), fftconvolve) instead of the v7 flat 0.5% on +-0.5 mm.
5. f-impurity (G3): fraction F_IMP = 5% split +-, screen shift DSCR_IMP = 1.97 um
   (conservative: includes pre-slit velocity accrual); impurity records counted noise-only.
6. m13: HONEST Stage-1 gate -- Sec.-5.4 grid stats over 24 points excluding x_c = 0
   (undefined relative error at the field zero; reported as absolute check) PLUS a new
   gravity-map identity gate at the actual Configuration-E geometry.
7. m16: seed-paired (CRN) uniform-vs-optimal allocation comparison inside stage5.
8. De-hazard: stage2_reference.py is now GATE-ONLY (v7 wrote colliding
   stage2_maps.pkl/stage2_results.pkl formats); record-noise guard fixed for float32.
9. CONFLICT-2 pairing ceiling: analysis-level, folded by collect_v8_summary.py into the
   campaign-tier durations via v8_config.eff_window(z) (trims top edge x3 at early planes;
   low edge unaffected).

## v8.1 changelog (post first local run, 2026-07-04)
Your run exposed one crash and a failing Gate 1; forensics found and fixed SIX defects
(each mechanism proven before the fix; details in STATE.md):
 f1. stage2_run NameError: the mixture patch had displaced the line defining the
     analyst-side model density. Restored (and the smoke suite now calls the actual
     functions rather than reimplementing them).
 f2. Cumulative left-boundary truncation: adaptive x_c grids started at the 1e-4 support
     edge, dropping N(x_c[0]) != 0. Integration domain now extends to a 1e-8 tail
     (GRID_REL) while the analysis validity mask stays at 1e-4.
 f3. x = 0 bin-edge convention mismatch: binsum is left-inclusive, ctrg was right-
     inclusive; the fine-grid point exactly at x = 0 (representable on both lattices)
     carried a one-bin lever error -- invisible in v7 (bin CENTER there; larger t),
     catastrophic against v8's 25x smaller deep-plane moment signal. side='right'.
 f4. O(h^2) mixed Simpson/trapezoid cumulative: h/lambda_f is nearly constant across
     planes under v8 kinematics, leaving a plane-independent ~1.8 um/s error. Replaced
     by cubic-spline antiderivative (sim_core.cum_integrate) -- an estimator upgrade at
     the paper's fixed 151-point design. Gate 1 now passes at 0.13% worst-plane.
 f5. Incoherent-component records: using the coherent a_w on rho_inc injected weak-value
     x envelope artifacts at fringe minima. Maps now build per-slit record products
     (a_w_inc); both consumers use the proper mixture mean.
 f6. Estimator normalization/self-consistency: data moments are divided by the known
     (1-F)(1-P_sc), and the subtraction/division density is the MEASURED window density
     V*rho_w + (1-V)*rho_w_inc (Sec.-7 DoG semantics), not the ideal coherent one.
Plus: float32 guard 1e-300 -> 1e-30 in stage3_maps (cosmetic warning), and stage3
maps pass a_w_inc / rho_w_inc through from the stage-2 reference maps (eta-independent
at (1-V)*eta order, same x_c grids) so stage 4/6 see the mixture keys.
TWO PHYSICS FINDINGS (P4-material, verify in this run):
 p1. Per-plane matched-filter SNR is genuinely ~x0.42 of the v7 numbers (plane 0,
     Na=1e4: 0.418 vs 1.144): the D3 w->w_eff correction shrinks rho_w (x0.71) and
     raises per-record noise normalization (x1.55). v7's SNR was inflated by the very
     convention error this campaign fixes. Expect 5-sigma atom thresholds ~x5.7 and
     Table-3 tier durations to move accordingly.
 p2. Deep planes carry intrinsically less moment signal under v8 kinematics
     (mom - rho_w ~ t, and t(z=57cm) = 8.1 ms); the estimator fixes above keep the
     identity gate at the 0.1% level regardless.

## How to run
    cd simcode_v8
    ./run_all_v8.sh          # ~20-40 min; everything logged to v8_run_log.txt
Python >= 3.10, numpy, scipy (matplotlib only for the fig scripts, not the chain).
All RNGs seeded; the run is bit-reproducible.

## What to send back
    v8_summary.json          (all headline numbers, machine-readable)
    v8_run_log.txt           (full console log incl. gates)
    sim/stage2_results.pkl, sim/stage3_bias_results.pkl, sim/stage4_results.pkl,
    sim/stage5_inject_Na1e6.pkl, sim/stage5_inject_Na1e7.pkl,
    sim/stage6_results.pkl, sim/fig4_data.pkl, sim/true_traj.pkl,
    sim/stage1_verification.pkl, sim/stage2ref_gate.pkl
KEEP LOCALLY (large, needed only if we iterate): sim/stage2_maps.pkl,
sim/stage3_maps_*.pkl, sim/stage4_design.pkl.

## Acceptance gates (chain aborts loudly if violated)
- v8_config self-checks pass (kinematics/couplings vs manuscript verification scripts)
- Stage-1a (Sec.-5.4 grid, 24 pts): max rel err < 4%   [smoke here: mean 1.04%, max 2.52%]
- Stage-1b (gravity map): max rel err < 4%             [smoke here: max 1.64%]
- Gate 1 identity (cumulative vs direct): worst plane < 5%
Note: make_fig2/3.py are v7-constants standalone scripts; they regenerate at P6, not here.
