"""Collect every manuscript-bound number from the v8 chain into v8_summary.json.
Sections keyed to the P4 rewrite manifest (PENDING_NUMBERS.md). Run from simcode_v8/."""
import numpy as np, pickle, json, sys
sys.path.insert(0, 'sim')
import v8_config as C

J = {}
zj = np.linspace(0.05, 0.57, 20)

def f(x): return float(x)
def arr(x): return [float(v) for v in np.asarray(x).ravel()]

# ---------- stage 1 gates ----------
s1 = pickle.load(open('sim/stage1_verification.pkl', 'rb'))
J['stage1'] = dict(grid24_mean_pct=f(s1['grid_mean']*100), grid24_max_pct=f(s1['grid_max']*100),
                   xc0_absmax_ms=f(s1['xc0_absmax']),
                   gravity_mean_pct=f(s1['gravity_mean']*100), gravity_max_pct=f(s1['gravity_max']*100),
                   pass_1a=s1['pass_1a'], pass_1b=s1['pass_1b'])
g1 = pickle.load(open('sim/stage2ref_gate.pkl', 'rb'))
devs = np.array([g1[k]['dev'] for k in range(20)])
J['gate1_identity'] = dict(worst_pct=f(devs.max()*100), median_pct=f(np.median(devs)*100))

# ---------- stage 2: raw pipeline ----------
s2 = pickle.load(open('sim/stage2_results.pkl', 'rb'))
rms = np.array([r['rms'] for r in s2['runs']]); env = np.array([r['env'] for r in s2['runs']])
J['stage2'] = dict(n_real=len(s2['runs']),
                   traj_rms_um=f(rms[:, -1].mean()*1e6), traj_rms_sd_um=f(rms[:, -1].std(ddof=1)*1e6),
                   traj_rms_pct_fringe=f(rms[:, -1].mean()/C.FRINGE*100),
                   clean_rms_um=f(s2['clean']['rms'][-1]*1e6),
                   clean_pct_fringe=f(s2['clean']['rms'][-1]/C.FRINGE*100),
                   env_ratio=f(env.mean()), env_ratio_sd=f(env.std(ddof=1)), env_clean=f(s2['clean']['env']))
xr = np.linspace(-950e-6, 950e-6, 381)
Vs = np.array([r['V'] for r in s2['runs']])
selr = np.abs(xr) < 400e-6
J['stage2']['sigv_z30_mms'] = f(Vs[:, 9, selr].std(axis=0).mean()*1e3)
J['stage2']['sigv_z5_mms']  = f(Vs[:, 0, selr].std(axis=0).mean()*1e3)

# ---------- stage 3: dipole bias ----------
s3 = pickle.load(open('sim/stage3_bias_results.pkl', 'rb'))
J['stage3'] = dict(zj_cm=arr(zj*100),
                   wiggle_scale_ums=arr([C.d/C.t_of_z(z)*1e6 for z in zj]),
                   wiggle_1em3_ums=arr(np.array(s3['wiggle_1em3'])*1e6),
                   D_1em3_ums=arr(np.array(s3['D_1em3'])*1e6),
                   D_1em2_ums=arr(np.array(s3['D_1em2'])*1e6),
                   sigv_1e4_ums=arr(np.array(s3['sigv_1e4'])*1e6))
wig = np.array([C.d/C.t_of_z(z) for z in zj])
for eb, er, Dv in [('1e-3', '1e-3', np.array(s3['D_1em3'])), ('1e-2', '1e-4', 0.1*np.array(s3['D_1em2'])),
                   ('1e-2', '1e-3', np.array(s3['D_1em2'])), ('1e-3', '1e-4', 0.1*np.array(s3['D_1em3']))]:
    fr = Dv/wig
    J['stage3'][f'cell_{eb}_{er}_worstfrac_pct'] = f(fr.max()*100)

# ---------- stage 4: GLS ----------
s4 = pickle.load(open('sim/stage4_results.pkl', 'rb'))
snr = np.array(s4['snr_1e4'])
J['stage4'] = dict(snr_1e4=arr(snr), snr_global=f(np.sqrt((snr**2).sum())),
                   Na_5sigma=arr((5.0/snr)**2*1e4))
J['stage4']['sweep'] = {}
for Na, dd in s4['sweep'].items():
    tr = np.array(dd['traj']); fr = np.array(dd['field_rel'])
    J['stage4']['sweep'][f'{Na:.0e}'] = dict(
        traj_um=f(tr.mean()*1e6), traj_sd_um=f(tr.std(ddof=1)*1e6),
        traj_pct_fringe=f(tr.mean()/C.FRINGE*100),
        field_rel_z5_pct=f(fr[:, 0].mean()*100), field_rel_z30_pct=f(fr[:, 9].mean()*100),
        field_rel_z57_pct=f(fr[:, 19].mean()*100))
J['stage4']['bias_gls_um'] = {f'{k[0]}_{k[1]}': f(v*1e6) for k, v in s4['bias_gls'].items()}
J['stage4']['bias_gls_pct_fringe'] = {f'{k[0]}_{k[1]}': f(v/C.FRINGE*100) for k, v in s4['bias_gls'].items()}

# ---------- stage 5: allocation ----------
# SE-smeared background in-window fraction per plane (from the stage-2 maps' bg_b); stored so that the master can
# check it from the shipped summary without the 60-MB maps (harness fix 2026-08-25). Recomputed whenever the chain runs.
try:
    _Mm = pickle.load(open('sim/stage2_maps.pkl', 'rb'))
    J['stage2']['se_inwindow_min_pct'] = f(min(float(np.asarray(m['bg_b']).sum()) for m in _Mm)*100)
    J['stage2']['se_inwindow_pct'] = arr([float(np.asarray(m['bg_b']).sum())*100 for m in _Mm])
except FileNotFoundError:
    pass
J['stage5'] = {}
for tag in ['1e6', '1e7']:
    try:
        s5 = pickle.load(open(f'sim/stage5_inject_Na{tag}.pkl', 'rb'))
    except FileNotFoundError:
        continue
    e = {'zeroth_gain': f(s5['zeroth'])} if 'zeroth' in s5 else {}
    if 'A' in s5:
        A = np.array(s5['A']); e['shares_opt'] = arr(np.sqrt(A)/np.sqrt(A).sum())
        e['additive_rms_um'] = f(np.sqrt(A.sum()/float(tag))*1e6)
    if 'A' in s5:   # additive-model optimal-allocation gain, exactly as stage5_alloc.py defines and prints it
        A = np.array(s5['A']); G = float(np.sqrt(20*A.sum())/np.sqrt(A).sum())
        e['additive_gain'] = f(G); e['additive_gain_pct'] = f(100*(G - 1))
    if 'm16' in s5:
        gn = np.array(s5['m16']['gain'])
        e['m16_gain_mean'] = f(gn.mean()); e['m16_gain_se'] = f(gn.std(ddof=1)/np.sqrt(gn.size))
        e['m16_benefit_pct'] = f((gn.mean()-1)*100)
    J['stage5'][tag] = e

# ---------- stage 6: field maps + campaign tiers (with CONFLICT-2 pairing ceiling) ----------
s6 = pickle.load(open('sim/stage6_results.pkl', 'rb'))
J['stage6'] = dict(osc={}, cells={})
for k, dd in s6.get('osc', {}).items():
    J['stage6']['osc'][str(k)] = dict(snr_osc_1e4=f(dd['snr_osc_1e4']),
                                      vscale_ums=f(dd['vscale']*1e6), osc_amp_ums=f(dd['osc_amp']*1e6))
for key, dd in s6.items():
    if isinstance(key, tuple):
        k, Na = key
        J['stage6']['cells'][f'p{k}_Na{Na:.0e}'] = dict(
            rel_pct=f(np.mean(dd['rel'])*100), rel_se_pct=f(np.std(dd['rel'], ddof=1)/np.sqrt(len(dd['rel']))*100),
            abs_ums=f(np.mean(dd['abs'])*1e6))
# Na(10%) per early plane by log-log fit of rel(Na); campaign durations vs eff_window
J['campaign'] = {}
for k in [0, 1, 2]:
    nas = sorted([key[1] for key in s6 if isinstance(key, tuple) and key[0] == k])
    rel = np.array([np.mean(s6[(k, Na)]['rel']) for Na in nas])
    lo, hi = C.eff_window(zj[k])
    fit = np.polyfit(np.log(nas), np.log(rel), 1)
    Na10 = float(np.exp((np.log(0.10) - fit[1])/fit[0]))
    Na15 = float(np.exp((np.log(0.15) - fit[1])/fit[0]))
    Na20 = float(np.exp((np.log(0.20) - fit[1])/fit[0]))
    J['campaign'][f'plane{k}'] = dict(
        z_cm=f(zj[k]*100), slope=f(fit[0]),
        Na_at_10pct=Na10, Na_at_15pct=Na15, Na_at_20pct=Na20,
        eff_window_per_day=[f(lo), f(hi)],
        ceiling_per_day=f(C.eff_window(zj[k])[1]),
        days_10pct=[f(151*Na10/hi), f(151*Na10/lo)],
        days_15pct=[f(151*Na15/hi), f(151*Na15/lo)])
J['campaign']['note'] = ('durations = 151 curtain positions x Na(target) / eff_window; '
                         'eff_window folds the CONFLICT-2 pairing ceiling (top edge only; '
                         'low edge unaffected, ceiling 9.6e6/day unselected > 3.2e6 source).')
J['config'] = dict(T_s=f(C.T), FRINGE_um=f(C.FRINGE*1e6), W_um=f(C.W*1e6),
                   SIG_ref=f(C.SIG_of(0.05)), SIG_last=f(C.SIG_of(0.57)),
                   V_ref=f(C.V_of(0.05)), Psc_ref=f(C.Psc_of(0.05)),
                   F_IMP=C.F_IMP, DSCR_IMP_um=f(C.DSCR_IMP*1e6), DU_PAIR_mms=f(C.DU_PAIR*1e3))

# ---------- fig4: GLS noiseless floor + sweep finals (guarded) ----------
try:
    F = pickle.load(open('sim/fig4_data.pkl', 'rb'))
    J['fig4'] = dict(clean_final_um=f(F['clean']['rms'][-1]*1e6),
                     clean_pct_fringe=f(F['clean']['rms'][-1]/C.FRINGE*100),
                     clean_env=f(F['clean']['env']),
                     TT_final_std_um=f(F['TT_final_std']*1e6))
    for Na in [1e4, 1e5, 1e6, 1e7]:
        if Na in F:
            rc = np.array(F[Na]['rms_curves'])
            J['fig4'][f'final_{Na:.0e}'] = dict(rms_um=f(rc[:, -1].mean()*1e6),
                                                sd_um=f(rc[:, -1].std(ddof=1)*1e6),
                                                env=f(np.mean(F[Na]['env'])))
except FileNotFoundError:
    J['fig4'] = 'fig4_data.pkl not present at collection time'

json.dump(J, open('v8_summary.json', 'w'), indent=1)
print("saved v8_summary.json")
print(json.dumps({k: (v if not isinstance(v, dict) else '...') for k, v in J.items()}, indent=1))
