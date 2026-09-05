"""P4 authoritative number set. Every value quoted in the P4 manuscript edits must appear here,
computed from the production pickles / v8_summary.json. Run from simcode_v8/."""
import numpy as np, pickle, json, sys
sys.path.insert(0,'sim'); import v8_config as C
J = json.load(open('v8_summary.json'))
print("== GATES ==")
print(f"S1a mean/max: {J['stage1']['grid24_mean_pct']:.2f}% / {J['stage1']['grid24_max_pct']:.2f}%  (24 pts, xc=0 excluded; |v(0)|<{J['stage1']['xc0_absmax_ms']:.1e} m/s)")
print(f"S1b gravity max: {J['stage1']['gravity_max_pct']:.2f}%   Gate-1 identity worst/median: {J['gate1_identity']['worst_pct']:.2f}% / {J['gate1_identity']['median_pct']:.2f}%")
print("== CLEAN FLOOR / TABLE 2 ==")
print(f"clean-pipeline final RMS: {J['stage2']['clean_rms_um']:.1f} um = {J['stage2']['clean_pct_fringe']:.2f}% of fringe; env ratio {J['stage2']['env_clean']:.3f}")
print("== DETECTION TIER (plane 0, 151 positions, pairing ceiling) ==")
snr0 = J['stage4']['snr_1e4'][0]
Na_pos = (5.0/snr0)**2*1e4; Ntot = 151*Na_pos
lo,hi = C.eff_window(0.05)
print(f"SNR(1e4, z=5cm) = {snr0:.4f} -> Na/position(5sig) = {Na_pos:.3e} -> TOTAL = {Ntot:.3e} atoms")
print(f"eff_window z=5cm = [{lo:.2e}, {hi:.2e}]/day -> DURATION = {Ntot/hi:.1f} - {Ntot/lo:.1f} days")
print(f"(vs v7 tier 2.87e7, 0.33-9.0 d: atom factor x{Ntot/2.87e7:.2f})")
print("== Na SWEEP (GLS pipeline; 1e4-1e5 traj non-informative/chaotic) ==")
for tag in ['1e+06','1e+07']:
    d=J['stage4']['sweep'][tag]
    print(f"Na={tag}: traj {d['traj_um']:.0f}+-{d['traj_sd_um']:.0f} um ({d['traj_pct_fringe']:.0f}% fringe); field z=5cm {d['field_rel_z5_pct']:.1f}%")
print("== DIPOLE BIAS (GLS end-to-end, % of fringe) + 5% crossing ==")
bg = J['stage4']['bias_gls_pct_fringe']
for k in ['1e-2_1e-3','1e-3_1e-3','1e-2_1e-4','1e-3_1e-4']:
    if k in bg: print(f"  eta,eps={k}: {bg[k]:.3f}%")
p6 = max(bg.get('1e-3_1e-3',0), bg.get('1e-2_1e-4',0)); p5 = bg.get('1e-2_1e-3',None)
if p5: 
    slope = np.log10(p5/p6)
    x5 = 1e-6*10**(np.log10(5.0/p6)/slope)
    print(f"  5%-of-fringe crossing (log-log, worst 1e-6 cell -> 1e-5 cell): product ~ {x5:.1e}  (v7-era ~2e-6; budget 1.4e-6 unchanged, more headroom)")
print("== ALLOCATION / m16 ==")
for tag in ['1e6','1e7']:
    e=J['stage5'][tag]
    print(f"Na={tag}: additive RMS {e.get('additive_rms_um',float('nan')):.0f} um; m16 CRN gain {e.get('m16_gain_mean',float('nan')):.3f} (+{e.get('m16_benefit_pct',float('nan')):.1f}% +- {100*e.get('m16_gain_se',float('nan')):.1f})")
print("== TRAJECTORIES: quote fig4 (20 realizations; heavy-tailed -- give mean+-SD) ==")
if 'fig4' in J and isinstance(J['fig4'], dict):
    F=J['fig4']
    print(f"GLS noiseless floor: {F['clean_final_um']:.2f} um ({F['clean_pct_fringe']:.2f}% fringe), env {F['clean_env']:.3f}")
    for tag in ['1e+06','1e+07']:
        d=F[f'final_{tag}']
        print(f"Na={tag}: {d['rms_um']:.0f} +- {d['sd_um']:.0f} um ({d['rms_um']/C.FRINGE/1e6*100:.0f}% fringe), env {d['env']:.2f}")
print("== m2: SE-smeared background in-window fraction per plane (from maps bg_b) ==")
M = pickle.load(open('sim/stage2_maps.pkl','rb')) if __import__('os').path.exists('sim/stage2_maps.pkl') else None
if M:
    fr=[m['bg_b'].sum() for m in M]
    print(f"  planes 0/9/19: {fr[0]:.3f} / {fr[9]:.3f} / {fr[19]:.3f}; min over planes {min(fr):.3f}")
else:
    print("  (rebuild sim/stage2_maps.pkl to compute)")
