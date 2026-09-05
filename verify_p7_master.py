"""NJP-120031 MASTER VERIFICATION (v8.1 chain; extended through the P8 referee-driven revision, 2026-08).
Run from the directory that contains this file (the supplementary-code root): it expects
  simcode_v8/ (chain + subprocess verifiers + v8_summary.json + sim/true_traj.pkl), work/main.tex,
  official_fig4_data.pkl (the official Fig-4 data anchor) and simcode/sim/stage4_results.pkl (a v7-era input
  read by one closure check). Every manuscript-quoted number is recomputed from closed forms / production
pickles / the official run and PASS/FAILed; later sections run structural sweeps of main.tex and the P8
subprocess verifiers (Appendix B, symbol table, stale-pattern/invariant sweep, review-fix pins).
ASCII output; exit 1 on any FAIL."""
import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))
import numpy as np, pickle, json, re, subprocess, sys

NP, NF = 0, 0
def ok(name, val, ref, tol=0.02, kind='rel'):
    global NP, NF
    good = (abs(val-ref) <= tol*abs(ref)) if kind=='rel' else (abs(val-ref) <= tol)
    NP += good; NF += (not good)
    print(f"{'PASS' if good else 'FAIL':4s} {name:64s} calc={val:12.5g}  quoted={ref:.5g}")
def okr(name, val, lo, hi):
    global NP, NF
    good = lo <= val <= hi
    NP += good; NF += (not good)
    print(f"{'PASS' if good else 'FAIL':4s} {name:64s} calc={val:12.5g}  band=[{lo:g},{hi:g}]")
def sec(t): print("\n" + "="*100 + f"\n{t}\n" + "="*100)

sys.path.insert(0, 'simcode_v8/sim')
import v8_config as C
TEX = open('work/main.tex').read()
J = json.load(open('simcode_v8/v8_summary.json'))
hbar=1.054571817e-34; h=6.62607015e-34; kB=1.380649e-23; g=9.80665; M=C.M

sec("A. v8 KINEMATICS / GEOMETRY (closed forms vs quoted)")
ok("fringe hT/(Md) (um)", h*C.T/(M*C.d)*1e6, 420.3, 0.001)
ok("total time T (s)", C.T, 0.229, 0.005)
ok("v_z(L) (m/s)", C.vz_of(0.60), 3.744, 0.001)
ok("release fall v0^2/2g (cm)", 1.5**2/(2*g)*100, 11.5, 0.01)
ok("kinetic temperature at v0 (mK)", 0.5*M*1.5**2/kB*1e3, 12, 0.05)
ok("coherence velocity hbar/(Md) (mm/s)", hbar/(M*C.d)*1e3, 0.29, 0.02)
ok("recoil velocity (mm/s)", h/(M*780.241e-9)*1e3, 5.9, 0.01)
# tab:kinematics five rows (z, vz, tau_a, t, d/t)
rows = [(0.05,1.797,1.113,0.0303,82.4),(0.0774,1.941,1.030,0.0450,55.6),
        (0.1047,2.075,0.964,0.0586,42.7),(0.2963,2.839,0.704,0.1366,18.3),
        (0.57,3.665,0.546,0.2207,11.3)]
for z, vzq, tq, ttq, wq in rows:
    ok(f"  tab:kin z={z*100:.2f}: v_z", C.vz_of(z), vzq, 0.002)
    ok(f"  tab:kin z={z*100:.2f}: tau_a (us)", C.tau_of(z)*1e6, tq, 0.005)
    ok(f"  tab:kin z={z*100:.2f}: t(z) (s)", C.t_of_z(z), ttq, 0.005)
    ok(f"  tab:kin z={z*100:.2f}: d/t wiggle (um/s)", C.d/C.t_of_z(z)*1e6, wq, 0.005)

sec("B. TABLE-2 COUPLINGS / STATISTICS (config anchors vs quoted)")
ok("g_a measuring (rad)", C.GA_M_REF, 3673, 0.001)
ok("g_a compensator (rad)", C.GA_C_REF, -3679, 0.001)
ok("Lambda_eff(z_ref)", -np.log(C.V_of(0.05)), 8.65e-3, 0.01)
ok("V(z_ref)", C.V_of(0.05), 0.991, 0.001)
ok("V^2 + D^2 with D=0.131", C.V_of(0.05)**2 + 0.131**2, 1.0, 0.002)
ok("P_sc(z_ref)", C.Psc_of(0.05), 0.0144, 0.01)
ok("SIG(z_ref) record noise", C.SIG_of(0.05), 6.146, 0.001)
ok("P8-8: per-shot SNR = 1/SIG(z_ref) quoted as 'of order 0.16' (Sec 10; response R1.20)", 1.0/C.SIG_of(0.05), 0.16, 0.02)
ok("tag significance sqrt(Na*Lam) @1e4", np.sqrt(1e4*8.65e-3), 9.3, 0.01)
ok("Na*(S_tag=10)", 100/8.65e-3, 1.16e4, 0.01)
ok("S_tag range low (last plane)", np.sqrt(1e4*8.65e-3*C.tau_of(0.57)/C.tau_of(0.05)), 7, 0.08)
snr0 = J['stage4']['snr_1e4'][0]
ok("Table-2 detection: Na/position (5sig, z_ref)", (5/snr0)**2*1e4, 1.4e6, 0.03)
ok("Table-2 detection: total (151 pos)", 151*(5/snr0)**2*1e4, 2.2e8, 0.03)
F = pickle.load(open('official_fig4_data.pkl','rb'))
tr7 = np.array(F[1e7]['rms_curves'])[:, -1]
ok("Table-2 sigma_traj @1e7 (mm) [OFFICIAL fig4]", tr7.mean()*1e3, 0.18, 0.03)
ok("  same, %% of fringe", tr7.mean()/420.3e-6*100, 44, 0.01)
# free-space single-pass reference (Table 3 reference block): Lambda_1 from the closed form of Eq. (freespace_ceiling)
# with the two-level, detuning-independent ratio phi_1^2/p_sc,1 = hbar omega Gamma/(4 pi w^2 I_sat)  [same constants as make_fig2.py]
_hbar = 1.054571817e-34; _om = 2*np.pi*384.230e12; _Gam = 2*np.pi*6.07e6; _Isat = 25.03; _w = 2e-6; _Psc_ref = 0.017
_ratio = _hbar*_om*_Gam/(4*np.pi*_w**2*_Isat)
ok("two-level phi_1^2/p_sc,1 (Fig 2a; det.-independent)", _ratio, 7.718e-3, 0.002)
LAM1 = (C.d/C.W)**2*np.exp(-C.d**2/(2*C.W**2))*_ratio*_Psc_ref
ok("free-space Lambda_1 [Table 3 ref block, P_sc = 0.017] (LIVE, harness fix 2026-08-25)", LAM1, 8.6e-5, 0.02)
ok("free-space Na(S_tag=10)", 100/LAM1, 1.2e6, 0.04)

sec("C. FIRST-LIGHT (power/30)")
ok("first-light g_a", C.GA_M_REF/30, 1.2e2, 0.03)
ok("first-light Lambda_eff", 8.65e-3/30, 2.9e-4, 0.02)
ok("first-light Na(S_tag=10)", 100/(8.65e-3/30), 3.5e5, 0.02)

sec("D. LEGACY VERIFY SCRIPTS (G1-G3 sets, run as subprocesses)")
# P8-4 (2026-08-24): verify_appB.py (Appendix B derivation chain; lives at the campaign root) added to the subprocess set
# P8-7g (2026-08-24): verify_glossary.py (Appendix C symbol table vs body) added to the subprocess set
# P8-8 (2026-08-24): verify_p8_8.py (stale-pattern sweep + do-not-regress invariants + bibliography cross-checks) added
# P8-8b (2026-08-25): verify_p8_8b.py (adversarial-review fix pins) added
for scr in ['verify_kinematics.py', 'verify_p2_couplings.py', 'verify_p3_systematics.py', 'verify_appB.py', 'verify_glossary.py', 'verify_p8_8.py', 'verify_p8_8b.py']:
    path = f'simcode_v8/{scr}' if scr not in ('verify_appB.py', 'verify_glossary.py', 'verify_p8_8.py', 'verify_p8_8b.py') else scr
    r = subprocess.run([sys.executable, path], capture_output=True, text=True)
    p = r.stdout.count('PASS'); f = r.stdout.count('FAIL')
    global_ok = (f == 0 and r.returncode in (0, None) and p > 0)
    NPF = 'PASS' if global_ok else 'FAIL'
    print(f"{NPF:4s} {scr:64s} sub-checks: {p} PASS / {f} FAIL")
    if not global_ok:
        NF += 1; print(r.stdout[-600:])
    else:
        NP += 1

sec("E. PRODUCTION / P4 NUMBER SET (summary + official log)")
ok("Gate 1a mean (%)", J['stage1']['grid24_mean_pct'], 1.04, 0.01)
ok("Gate 1a max (%)", J['stage1']['grid24_max_pct'], 2.52, 0.005)
ok("Gate 1b gravity max (%)", J['stage1']['gravity_max_pct'], 1.64, 0.01)
ok("Gate 1 identity worst (%)", J['gate1_identity']['worst_pct'], 0.16, 0.05)
ok("Gate 1 identity median (%)", J['gate1_identity']['median_pct'], 0.12, 0.05)
ok("clean-pipeline floor (um)", J['stage2']['clean_rms_um'], 1.9, 0.03)
ok("clean floor (% fringe)", J['stage2']['clean_pct_fringe'], 0.45, 0.03)
okr("direct pipeline @1e4 (mm) [chaotic cell; shipped run of record 12.78; prose '~13'] (LIVE, harness fix 2026-08-25)", J['stage2']['traj_rms_um']/1000, 10.0, 17.0)
ok("direct pipeline @1e4 rounds to the quoted ~13 mm", round(J['stage2']['traj_rms_um']/1000), 13, 0.0)
snr = np.array(J['stage4']['snr_1e4'])
ok("SNR ladder z=5cm", snr[0], 0.42, 0.01)
ok("SNR ladder z=13cm", snr[3], 0.06, 0.09)
ok("global pooled SNR", J['stage4']['snr_global'], 0.47, 0.01)
ok("front-load: planes 0-2 share (%)", np.sqrt((snr[:3]**2).sum())/J['stage4']['snr_global']*100, 98, 0.005)
for k, q in [(0, 1.4e6), (1, 8.2e6), (2, 2.8e7)]:
    ok(f"5sig threshold plane {k}", J['stage4']['Na_5sigma'][k], q, 0.03)
lo, hi = C.eff_window(0.05)
Ntot = 151*(5/snr0)**2*1e4
ok("detection duration low (d)", Ntot/hi, 7.5, 0.02)
ok("detection duration high (d)", Ntot/lo, 67, 0.01)
ok("flux window width (x)", hi/lo, 9, 0.01)
bg = J['stage4']['bias_gls_pct_fringe']
ok("dipole 1e-7 cell (% fringe)", bg['1e-3_1e-4'], 0.214, 0.02)
ok("dipole 1e-6 cell A", bg['1e-3_1e-3'], 1.945, 0.02)
ok("dipole 1e-6 cell B", bg['1e-2_1e-4'], 1.164, 0.02)
ok("dipole 1e-5 cell", bg['1e-2_1e-3'], 8.939, 0.02)
x5 = 1e-6*10**(np.log10(5.0/max(bg['1e-3_1e-3'],bg['1e-2_1e-4']))/np.log10(bg['1e-2_1e-3']/max(bg['1e-3_1e-3'],bg['1e-2_1e-4'])))
ok("5%%-crossing product", x5, 4e-6, 0.10)
ok("cell ratio (factor 1.7)", bg['1e-3_1e-3']/bg['1e-2_1e-4'], 1.7, 0.02)
cN10 = J['campaign']['plane0']['Na_at_15pct']*151
ok("15%% map total atoms (plane 0)", cN10, 4e9, 0.05)
ok("15%% map days at ceiling", J['campaign']['plane0']['Na_at_15pct']*151/hi, 130, 0.05)
ok("plane-1/plane-0 cost (x5)", J['campaign']['plane1']['Na_at_15pct']/J['campaign']['plane0']['Na_at_15pct'], 5, 0.05)
ok("plane-2/plane-0 cost (x16)", J['campaign']['plane2']['Na_at_15pct']/J['campaign']['plane0']['Na_at_15pct'], 16, 0.05)
osc = J['stage6']['osc']
ok("osc amp high (um/s)", osc['0']['osc_amp_ums'], 14, 0.05)
ok("osc amp low (um/s)", osc['2']['osc_amp_ums'], 10, 0.02)
ok("field scale (~360 um/s)", osc['0']['vscale_ums'], 360, 0.02)
ok("osc 5sig threshold (>=2e9)", (5/osc['0']['snr_osc_1e4'])**2*1e4, 2e9, 0.25)
fg = J['fig4']
ok("fig4 1e4 (mm)", fg['final_1e+04']['rms_um']/1e3, 1.25, 0.01)
ok("fig4 1e4 sd (mm)", fg['final_1e+04']['sd_um']/1e3, 0.84, 0.01)
ok("fig4 1e5 (mm)", fg['final_1e+05']['rms_um']/1e3, 0.64, 0.01)
ok("fig4 1e6 (mm)", fg['final_1e+06']['rms_um']/1e3, 0.51, 0.01)
ok("fig4 1e7 (um)", fg['final_1e+07']['rms_um'], 184, 0.005)
ok("fig4 1e7 sd (um)", fg['final_1e+07']['sd_um'], 93, 0.01)
ok("fig4 1e7 %% fringe", fg['final_1e+07']['rms_um']/420.3, 0.44, 0.01)
ok("fig4 env @1e4", fg['final_1e+04']['env'], 7.3, 0.01)
ok("fig4 env @1e7", fg['final_1e+07']['env'], 1.5, 0.02)
ok("GLS noiseless floor (um)", fg['clean_final_um'], 3.4, 0.02)
ok("GLS floor (%% fringe)", fg['clean_pct_fringe'], 0.8, 0.01)
ok("GLS floor env", fg['clean_env'], 1.03, 0.005)
m16 = J['stage5']
okr("m16 @1e6 [shipped run of record 1.508; quoted 1.51+-0.31] (LIVE, harness fix 2026-08-25)", m16['1e6']['m16_gain_mean'], 1.45, 1.55)
ok("m16 @1e6 rounds to the quoted 1.51 (se 0.31)", round(m16['1e6']['m16_gain_mean'], 2), 1.51, 0.0)
ok("m16 @1e6 se rounds to the quoted 0.31", round(m16['1e6']['m16_gain_se'], 2), 0.31, 0.0)
ok("m16 @1e7", m16['1e7']['m16_gain_mean'], 1.11, 0.01)
ok("alloc additive gain @1e6 (G = 1.076 -> 7.6%) (LIVE, harness fix 2026-08-25)", m16['1e6']['additive_gain_pct'], 7.6, 0.02)
okr("additive-model gains at 1e6 and 1e7 inside the quoted 5-8%", min(m16['1e6']['additive_gain_pct'], m16['1e7']['additive_gain_pct']), 5.0, 8.0)
okr("  (upper)", max(m16['1e6']['additive_gain_pct'], m16['1e7']['additive_gain_pct']), 5.0, 8.0)
ok("SE in-window fraction min (>=99%%) [from the shipped summary; recomputed from the maps by collect_v8_summary.py]", J['stage2']['se_inwindow_min_pct'], 99.1, 0.005)

sec("F. FIGURE-2/3 DERIVED ANCHORS")
Gam=2*np.pi*6.07e6; Isat=25.03; w=2e-6; om=2*np.pi*3.8423e14
ok("fig2 ratio phi1^2/psc1", hbar*om*Gam/(4*np.pi*w**2*Isat), 7.7e-3, 0.01)
ok("fig2 ideal bound sig_res/2w^2", 3*(780e-9)**2/(2*np.pi)/(2*w**2), 0.036, 0.01)
ok("fig3 star SNR", np.sqrt(1e4*8.65e-3), 9.3, 0.005)
ok("fig3 single-pass Na", 100/LAM1, 1.2e6, 0.04)

sec("G. RECORDED-NOT-QUOTED (removed claims, best-effort rederivations)")
Pin = C.NG_M_REF*hbar*om/C.tau_of(0.05)   # input power per measuring arm (Ng is per-transit at input; F multiplies coupling)
print(f"REC  per-arm input power N_g*hbar*om/(F*tau_a) (measuring): {Pin*1e3:.2f} mW  [caption claim removed; record only]")
import os as _os2
if not _os2.path.exists('simcode_v8/sim/stage2_maps.pkl'):
    print('REC  stationary-phase residue: needs simcode_v8/sim/stage2_maps.pkl (rebuild the chain); record-only, no check affected')
else:
  Mm = pickle.load(open('simcode_v8/sim/stage2_maps.pkl','rb'))
  mp = Mm[0]; xf = 0.5*(np.linspace(-1.25e-3,1.25e-3,2501)[1:]+np.linspace(-1.25e-3,1.25e-3,2501)[:-1])
  prodb = mp['a_w'].astype(float)*mp['rho_b'][None,:]
  lev = xf[None,:]-mp['xc'][:,None]
  mom = ((lev*prodb).sum(axis=1)+mp['D_corr'])/C.W
  resid = np.abs(mom - mp['rho_w']); gross = (np.abs(lev*prodb)).sum(axis=1)/C.W
  core = mp['rho_w'] >= 0.05*mp['rho_w'].max()
  print(f"REC  stationary-phase residue |mom-rho_w|/gross (plane-0 core median): {np.median(resid[core]/gross[core])*100:.2f}%  [0.2% claim removed]")

sec("H. ADVERSARIAL-FINDINGS CLOSURE (F1-F12 grep evidence in main.tex)")
def gre(name, pattern, want_present=True, cnt=None):
    global NP, NF
    n = len(re.findall(pattern, TEX))
    good = (n == cnt) if cnt is not None else ((n > 0) == want_present)
    NP += good; NF += (not good)
    print(f"{'PASS' if good else 'FAIL':4s} {name:64s} matches={n}" + (f" (want {cnt})" if cnt is not None else ""))
gre("F1 gravity: fringe 420 present", r'420')
gre("F1 gravity: '735' only as historic/none", r'735', cnt=0)
gre("F1 gravity: v_z ~ 1 m/s relic absent", r'v_z \\sim 1 m/s', False)
gre("F2 Coriolis section present", r'Coriolis')
gre("F3 convention: w_eff defined", r'w_\\text\{eff\}')
gre("F4 sign convention note present", r'sign', True)
gre("F6 Ramos cited", r'Ramos2020')
gre("F6 Chapman cited", r'Chapman1995')
gre("F6 Aljunid cited", r'Aljunid2009')
gre("F7 Madelung corrected title", r'hydrodynamischer Form')
gre("F7 chimeric Madelung title absent", r'Maxwell-Lorentz-Einstein', False)
gre("F8 462uW appears only as the single-pass reference identifier (ref-block, Fig-2 caption, Table-3 caption) [amended P8-7 per ratified A3]", r'462', cnt=3)
gre("F9 pm0.5-mm window absent", r'\\pm 0\.5\}?-mm', False)
gre("F10 tier cost relic absent", r'2\.6\\times.*4\.2\\times', False)
gre("F11 crossing at 4e-6", r'crosses the \$5\\%\$ level near \$4\\times 10\^\{-6\}\$')
gre("F12 Eq17 bound 0.036 present", r'0\.036')
gre("F12 0.073 absent", r'0\.073', False)

sec("I. STRUCTURAL")
cites = set(re.findall(r'\\cite\{([^}]*)\}', TEX))
cites = set(c.strip() for grp in cites for c in grp.split(','))
bibs = set(re.findall(r'\\bibitem\{([^}]*)\}', TEX))
miss = cites - bibs; orph = bibs - cites
gre("all \\cite have \\bibitem", 'x', cnt=None) if False else None
good = (len(miss) == 0); NP += good; NF += (not good)
print(f"{'PASS' if good else 'FAIL':4s} {'all cite keys have bibitems':64s} missing={sorted(miss)}")
good = (len(orph) == 0); NP += good; NF += (not good)
print(f"{'PASS' if good else 'FAIL':4s} {'no orphan bibitems':64s} orphans={sorted(orph)}")
for pat, nm in [(r'tab:tiers', 'no tab:tiers refs'), (r'Tier[- ]', 'no Tier nomenclature'),
                (r'4\.5\\times 10\^\{-7\}', 'no stale budget'), (r'4\.4\\times 10\^3', 'no stale g_a'),
                (r'0\.989', 'no stale V'), (r'28-point', 'no 28-point claims'),
                (r'staged programme', 'no staged-programme phrase'),
                # P8-2 (2026-08-24): notation-collision bans, extended per the ratified P8-1 decision 2
                (r'\\Delta\(x_c', 'P8-2: no Delta(x_c,z_j) moment (renamed script-A)'),
                (r'\\widehat\{D\}', 'P8-2: no hatted-D denominator (renamed rho_w-hat)'),
                (r'N/\(t \\cdot D\)', 'P8-2: no N/(t.D) denominator form'),
                (r'V\(\\hat x\)', 'P8-2: no bare-V potential in App-A (renamed V_dipole)'),
                (r'\(§[0-9]', 'P8-2: no hard-coded section numbers')]:
    gre(f"sweep: {nm}", pat, False)
gre("P8-2: script-A moment present in Eq. numerator", r'\\mathcal\{A\}\(x_c, z_j\)')
# P8-3 (2026-08-24): the one new manuscript number of the phase -- resonant recoil vs fringe momentum, Sec. 1
ok("P8-3: recoil ratio d/lambda_D2 quoted as 3.2 h/d (Sec 1)", C.d/C.lamD2, 3.2, 0.01)
gre("P8-3: '3.2\\,h/d' present in Sec 1", r'3\.2\\,h/d')
gre("P8-3: old paragraph heading gone", r'Photon transverse momentum kick', False)
# P8-5 (2026-08-24): decoherence section rigor (R2-9)
gre("P8-5: explicit partial trace Tr_gamma present in Sec 9", r'\\mathrm\{Tr\}_\\gamma')
gre("P8-5: eq:partial_trace labelled once", r'\\label\{eq:partial_trace\}')
gre("P8-5: no chi_{L,R} conditional-state notation left", r'\\chi_', False)
gre("P8-5: no script-D Englert distinguishability in Sec 1 (plain D reserved for Englert)", r'\\mathcal\{D\}\^2', False)
gre("P8-5: 'overlap exponent' clarification present", r'overlap exponent')
# P8-6 (2026-08-24): every number quoted in the Fig-4 commentary / captions, from the OFFICIAL fig4 pickle + rebuilt true_traj
_F6 = pickle.load(open('official_fig4_data.pkl', 'rb')); _TT6 = pickle.load(open('simcode_v8/sim/true_traj.pkl', 'rb'))
_fin = np.array(_F6[1e7]['rms_curves'])[:, -1]*1e6; _z6 = _F6['zeval']*100
ok("P8-6: 1e7 final-RMS range min (um), quoted 82", _fin.min(), 82.4, 0.01)
ok("P8-6: 1e7 final-RMS range max (um), quoted 353", _fin.max(), 353.0, 0.01)
ok("P8-6: 1e7 final-RMS median (um), quoted 149", float(np.median(_fin)), 149.1, 0.01)
ok("P8-6: 1e7 mean RMS / noiseless floor > 50 ('more than fifty times')", _fin.mean()/(_F6['clean']['rms'][-1]*1e6), 54.5, 0.02)
_Vr, _Tr = _F6['rep_1e7']; _Tr = _Tr.astype(float)
_rms_r = np.sqrt(np.mean((_Tr - _TT6)**2, axis=0))*1e6; _disp = np.mean(_Tr - _TT6, axis=0)*1e6; _env = _TT6.std(axis=0)*1e6
_i20 = int(np.argmin(np.abs(_z6 - 20)))
ok("P8-6: representative 1e7 realization RMS at z = 20 cm (um), quoted 21", _rms_r[_i20], 21.0, 0.03)
ok("P8-6: true fan width (std) at z = 20 cm (um), quoted 45", _env[_i20], 45.3, 0.02)
ok("P8-6: representative realization: max mean displacement of the bundle (um), quoted 'up to ~300'", _disp.max(), 288.2, 0.02)
ok("P8-6: fraction of trajectories ending beyond +300 um, quoted 'about a quarter'", float((_Tr[:, -1] > 300e-6).mean()), 0.27, 0.02)
ok("P8-6: representative realization: max excursion (um), quoted '~ +500'", _Tr.max()*1e6, 530.0, 0.02)
ok("P8-6: representative realization = r0 with final RMS close to the 20-realization mean (um)", _rms_r[-1], 187.2, 0.01)
gre("P8-6: Fig-1 caption carries the dispersive/elastic clause", r'dispersive and elastic')
gre("P8-6: Fig-5 caption states same configuration as Fig-4", r'same configuration and pipeline as Figure')
gre("P8-6: 'Reading Figure' commentary paragraph present", r'\\paragraph\{Reading Figure')
# P8-7 (2026-08-24): abstract & framing (R1-9, R1-7, R2-10, A1, A3)
gre("P8-7: no configuration letters in prose (R1-7)", r'Configuration[~ -][A-F]\\b', False)
gre("P8-7: the one surviving letter: 'design-study point E' exactly once", r'design-study point E', True, 1)
gre("P8-7: no 'Config.-D' abbreviation", r'Config\\.-D', False)
gre("P8-7: A1 wording present (contributions + Sec 11)", r'verified numerically against the exact \(closed-form\) Bohmian field', True, 2)
gre("P8-7: no 'analytically verified' / 'verified analytically'", r'analytically verified|verified analytically', False)
gre("P8-7: Discussion mirror of contributions present", r'Set against the six contributions')
gre("P8-7g: Appendix C label present once", r'\\label\{app:notation\}', True, 1)
# P8-8 review F2 (2026-08-25): the fringe spacing is quoted uniformly as 420 um (3 s.f. of hT/(Md) = 420.3); logged in STATE
gre("P8-8/F2: Eq. (2) displays the fringe as 420 um", r'= 420~\\mu\\text\{m\}')
gre("P8-8/F2: no '420.3' anywhere in the manuscript (uniform 3-s.f. rounding)", r'420\.3', False)
# P8-8 review F1 (2026-08-25): float placement -- Fig 1 on a float page, Figs 3-5 with [!t]
gre("P8-8/F1: Fig 1 uses a float page", r'\\begin\{figure\}\[p\]\s*\\centering\s*\\includegraphics\[width=0\.92\\textwidth\]\{photon_curtain_fig1', True, 1)
_ab = TEX[TEX.find(r'\begin{abstract}')+16:TEX.find(r'\end{abstract}')]
_abw = max(len(re.sub(r'\$[^$]*\$', 'X', _ab.replace(r'\noindent', '')).split()), len(_ab.replace(r'\noindent', '').split()))   # P8-8 F3: the STRICTER of the two counting conventions
ok("P8-7: abstract word count <= 300 (IOP/NJP rule; counted with each math token as one word)", float(_abw <= 300), 1.0, 0.0)
print(f"      abstract word count = {_abw}")
ok("P8-7: abstract contains no \\ref/\\eqref/\\cite (IOP/NJP rule)", float(re.search(r'\\(ref|eqref|cite)\b', _ab) is None), 1.0, 0.0)
gre("P8-7: abstract free of internal notation (Lambda_eff, eta_bal, varphi_1, Gamma/Delta)", r'begin\{abstract\}[^\x00]*?(Lambda_\\text\{eff\}|eta_\\text\{bal\}|varphi_1|Gamma/\\Delta)[^\x00]*?end\{abstract\}', False)

print("\n" + "="*100)
print(f"MASTER VERIFICATION: {NP} PASS / {NF} FAIL")
sys.exit(0 if NF == 0 else 1)
