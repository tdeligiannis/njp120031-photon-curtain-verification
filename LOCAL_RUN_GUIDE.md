# Running and understanding the NJP-120031 scripts locally
(Linux terminal · Spyder · Claude Code) — written for the `supplementary_code/` folder of
`NJP120031_P8_overleaf_bundle_with_code.zip`, 2026-08-25.

## 0. The one-paragraph picture

Everything in the paper's numbers comes from one chain of Python scripts. `v8_config.py` holds the physical
constants, the free-fall kinematics and the per-plane couplings (the "operating configuration"; called
`Configuration E` in the code). `sim_core.py` knows how to build the exact double-slit wavefunction, its exact
guidance (Bohmian) velocity field, and the noiseless version of the paper's reconstruction formula, Eq. (18).
The `stage2…stage6` scripts are the forward model of the *experiment*: they simulate what the deflection
records would look like at each of the 20 curtain planes, add the record noise, run the estimator on the
noisy records, integrate trajectories through the reconstructed field, and measure how far they are from the
true ones. `make_fig*.py` draw Figures 2–5 from those results. The `verify_*.py` scripts are the independent
check: they recompute every number that the manuscript quotes and compare it to the text — the master
script, `verify_p7_master.py`, prints one PASS/FAIL line per number (178 at the moment) and exits non-zero
if anything fails. The only
large "data" shipped is `official_fig4_data.pkl`, the trajectory statistics of the authors' official run,
which the manuscript quotes and which the chain can regenerate to within the few-percent spread of a chaotic
regime; the small result files are the run of record for the two other chaotic-regime numbers.

## 1. Set up once

```bash
# 1. unpack somewhere with ~1 GB free (the chain writes ~740 MB of intermediates)
unzip NJP120031_P8_overleaf_bundle_with_code.zip
cd NJP120031_P8_overleaf_bundle_with_code/supplementary_code

# 2. a clean environment (conda shown; venv works the same way)
conda create -n njp python=3.12 numpy scipy matplotlib sympy mpmath -y
conda activate njp
python -c "import numpy, scipy, sympy, mpmath, matplotlib; print(numpy.__version__, scipy.__version__)"
```
Reference versions of the shipped results: Python 3.12.3, numpy 2.4.4, scipy 1.17.1, matplotlib 3.10.8,
sympy 1.14, mpmath 1.3. Newer versions are fine; all deterministic gates should reproduce to the digits
shown below regardless of machine. RAM: < 3 GB. Disk: ~740 MB after the full chain.

Working directory matters. Two rules cover everything:
- the four `verify_*.py` at the root and the master locate their inputs relative to their own file, so they
  run from anywhere; `verify_p2_couplings.py` must be run **from the root**
  (`python3 simcode_v8/verify_p2_couplings.py`) because it reads `simcode/sim/stage4_results.pkl`;
- every chain script, every figure script and `verify_p4_numbers.py` must be run **from `simcode_v8/`**
  (they open `sim/...` paths and import `sim/v8_config.py`).

## 2. Three ways to use the package

### Mode A — verify the manuscript against the shipped results (2 minutes)
```bash
cd supplementary_code
python3 verify_glossary.py          # Appendix C table vs body           -> "4 checks passed, 0 failed"
python3 verify_p8_8.py              # stale phrases / invariants / bib   -> "51 checks passed, 0 failed"
python3 verify_p8_8b.py             # numbers of the final revision pass -> "24 checks passed, 0 failed" (~1 min)
python3 verify_p7_master.py         # the full master on the shipped files -> "178 PASS / 0 FAIL" (~10 min)
cd simcode_v8
python3 verify_p4_numbers.py        # the Secs. 10-11 number set from the shipped result files
python3 make_fig2.py && python3 make_fig3.py && python3 make_fig5_new.py    # Figs 2, 3, 5 (PDF + *_check.png)
```
What you should see at the end of `make_fig2.py`: `ratio=7.718e-03  ideal=0.0363 ... crossover=0.35 GHz`; of
`make_fig3.py`: `SNR(star) = 9.30`. The regenerated PDFs reproduce the pinned numbers on any stack and are pixel-identical to the ones in the
Overleaf root under the reference matplotlib 3.10.8 (newer versions change bytes, not content).

### Mode B — regenerate everything and run the full verification (10–40 minutes)
```bash
cd supplementary_code
bash run_chain_and_verify.sh 2>&1 | tee my_reproduction.log
```
This deletes `simcode_v8/sim/*.pkl` (the chain refuses to start over stale pickles), runs
`simcode_v8/run_all_v8.sh` stage by stage, re-installs the official Fig-4 anchor, re-collects the summary
and runs the master. The last line must be `MASTER VERIFICATION: 178 PASS / 0 FAIL`.

If you prefer to watch each stage (or your terminal kills background jobs), run the stages yourself from
`simcode_v8/` in this order — each one writes its own pickle and the later ones read it:
```bash
python3 sim/v8_config.py                    # self-checks of constants, kinematics, couplings  (seconds)
python3 sim/sim_core.py                     # Stage-1 gates: 24-point grid + gravity map        (seconds)
python3 sim/stage2_maps.py 0 20             # record maps at the 20 planes -> stage2_maps.pkl  (~1 min)
python3 sim/stage2_reference.py             # Gate 1: estimator identity, plane by plane       (seconds)
python3 sim/stage2_run.py 30                # 30 direct-pipeline noise realizations             (~1 min)
for m in 0 1 2 3 4; do python3 sim/stage3_maps.py $m 0 20; done   # dipole-residual masks    (~2 min each)
python3 sim/stage3_analysis.py              # bias cells, budget test                          (seconds)
python3 sim/stage4_gls.py                   # GLS design matrices + noiseless gate             (~1 min)
python3 sim/stage4_run.py                   # detection SNRs, dipole bias through GLS          (~2 min)
python3 sim/stage5_alloc.py 1e6 12          # allocation study, Na = 1e6 (12 realizations)     (~1 min)
python3 sim/stage5_alloc.py 1e7 12          # ... and 1e7                                       (~1 min)
python3 sim/stage6_fieldmap.py              # early-plane field maps, campaign tiers           (seconds)
python3 sim/make_fig4_data.py               # trajectory statistics for Figs 4-5 (regenerated!) (~2 min)
cp ../official_fig4_data.pkl sim/fig4_data.pkl   # put the official anchor back (see Sec. 4 of this guide)
python3 collect_v8_summary.py               # -> v8_summary.json
python3 make_fig4_new.py                    # Fig 4 now possible (needs stage4_design.pkl from stage4_gls)
cd .. && python3 verify_p7_master.py        # the 178 checks
```
Numbers to expect from the deterministic gates (they reproduce exactly, on any machine):
`gravity-gate: mean 0.73% max 1.64%`, `STAGE-1 GATES: 1a PASS 1b PASS` (grid: mean 1.04 %, max 2.52 %),
`GATE 1 worst plane: 0.16% median plane: 0.12%`, `clean-pipeline limit: 1.9 um`, GLS gate `0.39 / 0.12 /
0.10 / 0.11 / 1.04 %`, `global (all planes pooled): 0.47`, `m16 ... 1.113 +- 0.192` at 1e7.
Numbers that legitimately vary by a few percent between machines (chaotic regime): the direct-pipeline
RMS at 1e4 (shipped run of record 12.78 mm, quoted as ~13 mm; the official machine's log shows 13.8 mm for a
different set of 30 realizations), the allocation ratio at 1e6 (shipped 1.508, quoted 1.51 +- 0.31), and the
`make_fig4_data.py` statistics (anchor 1252.8 um at 1e4 and 183.8 um at 1e7; a test machine gave 1237 and
183.9). The master reads the first two live from the shipped summary with chaos-honest bands and pins the
Fig-4 statistics to the anchor, which is why the anchor is re-installed.

### Mode C — check one thing you care about
- "Is Appendix B right?" → `python3 verify_appB.py` (~10 min; 35 checks: sympy identities + numerical
  propagation of the joint atom–meter state).
- "Where does number X in the paper come from?" → `grep -n "X" verification_log_master.txt` shows the
  PASS line with `calc=` (recomputed) and `quoted=` (manuscript); the line's name says which formula or
  pickle produced it; the same name appears in `verify_p7_master.py`.
- "What if I change a design parameter?" → edit `simcode_v8/sim/v8_config.py`, run `python3 sim/v8_config.py`
  (its self-checks will now FAIL against the manuscript's numbers, which is expected), then rerun the chain.
  The verifiers compare against the *manuscript*; a design change is an exploration, not a verification.

## 3. Environment-specific notes

### Linux terminal
Everything above is terminal-native. Use `tee` to keep logs (`... 2>&1 | tee stage.log`). Long stages run
in the foreground; if you must background them, use `nohup ... &` in a terminal that stays open, or `tmux`/
`screen`. Exit codes: every verifier exits 0 on all-pass, 1 otherwise, so they can be chained with `&&`.

### Spyder
1. Open Spyder in the `njp` environment (or set the interpreter under Preferences → Python interpreter).
2. Set the working directory: in the IPython console, `%cd /path/to/supplementary_code/simcode_v8` for chain
   and figure scripts, `%cd /path/to/supplementary_code` for the root verifiers (Spyder's Run → Configuration
   per file → "working directory" does the same thing per script).
3. Scripts that take command-line arguments (`stage2_maps.py 0 20`, `stage3_maps.py m 0 20`,
   `stage5_alloc.py 1e6 12`, `stage2_run.py 30`) must receive them: either Run → Configuration per file →
   "Command line options", or in the console `runfile('sim/stage3_maps.py', args='0 0 20', wdir='.')`, or
   simply use the console as a shell: `!python3 sim/stage3_maps.py 0 0 20`. Running such a script with F5
   and no arguments raises `IndexError` on `sys.argv` — that is the symptom.
4. The figure scripts call `matplotlib.use('Agg')` and save files; they will not pop up figures in the plot
   pane. Open the saved `photon_curtain_fig*.pdf` or the `fig*_check.png` next to it.
5. `verify_p7_master.py` runs its subprocess verifiers with `sys.executable`, so it uses Spyder's
   interpreter; run it with F5 (it changes into its own directory) and read the console.
6. The Variable Explorer is useful for the pickles: `import pickle; F = pickle.load(open('sim/fig4_data.pkl','rb'))`
   gives a dict keyed by the atom numbers (`1e4 … 1e7`), each holding `rms_curves` (20 realizations × 60 depths,
   metres), the `env` ratios, and `zeval`; `'clean'` is the noiseless floor, `'rep_1e7'` the representative
   realization drawn in Fig. 4(b).

### Claude Code (or any agentic terminal)
Give it the folder and the README and ask for one of the modes, e.g.
"Read supplementary_code/README.md, create a conda environment as it specifies, run Mode A and report the last
line of each verifier; then run run_chain_and_verify.sh in the foreground and report every gate line and the
master's final line." Two things to insist on: run long stages in the foreground (background jobs die when a
shell exits), and do not let it "fix" a FAIL by editing a verifier — a FAIL is information; the diagnosis is
whether the recomputed value (`calc=`) or the manuscript (`quoted=`) is wrong. Ask it to `grep FAIL` the log
and to compare regenerated pickles with the shipped ones (`sim/*.pkl`) before drawing conclusions.

## 4. Script atlas (what each file does; read this with the data-flow line)

Data flow:
`v8_config` → `sim_core` → `stage2_maps` → (`stage2_reference` gate, `stage2_run`) → `stage3_maps` ×5 →
`stage3_analysis` → `stage4_gls` → `stage4_run` → `stage5_alloc` → `stage6_fieldmap` → `make_fig4_data` →
`collect_v8_summary` → figures / verifiers.

| file | role | reads → writes |
|---|---|---|
| `sim/v8_config.py` | Single source of truth: constants (⁸⁷Rb, D1/D2), the vertical free-fall kinematics `vz(z)`, `t_of_z(z)`, `T`, the 20 plane depths, the two curtain arms' couplings anchored at z_ref = 5 cm and their per-plane scaling ∝ τ_a(z), the record-noise scale `SIG_of(z)`, the effective source-flux window. Running it executes its self-checks against the manuscript's verification scripts. | — → prints |
| `sim/sim_core.py` | `make_psi(d, σ₀, tmap)` returns the exact two-Gaussian wavefunction ψ(x,z), its derivative and the exact guidance field v(x,z); `recon_noiseless` evaluates Eq. (18) numerically (kernel width w_eff, plug-in density); `propagate`, `cum_integrate`. Running it executes the Stage-1 gates: the 24-point Sec.-5.4 grid and the gravity-map identity. | — → `stage1_verification.pkl` |
| `sim/stage2_maps.py lo hi` | For planes lo…hi builds the *record maps*: on an adaptive 151-point x_c grid, the weak-value profile a_w(x_f; x_c), the curtain-blurred density ρ_w and the per-slit components needed for the visibility mixture — i.e. what a noiseless experiment would record. | `v8_config`, `sim_core` → `stage2_maps.pkl` (60 MB) |
| `sim/stage2_reference.py` | Gate 1: the cumulative-estimator path (Eq. 29) against the direct formula (Eq. 18), plane by plane; must agree to < 5 % (it agrees to 0.16 %). | `stage2_maps.pkl` → `stage2ref_gate.pkl` |
| `sim/stage2_run.py n` | n noise realizations of the *direct* pipeline at 10⁴ atoms per point at the reference cell: adds record noise, reconstructs the field, integrates 200 trajectories, measures the RMS deviation from the true trajectories (which it also writes). Checkpoints every 5 realizations. | maps → `stage2_results.pkl`, `true_traj.pkl` |
| `sim/stage3_maps.py m lo hi` | The dipole back-action study: record maps under a residual dipole phase mask η ∈ {0, 10⁻³, 10⁻²} and two ε-offset masks (index m = 0…4) on a fine ±6 mm grid. | maps → `stage3_maps_eta*.pkl` (92 MB each) |
| `sim/stage3_analysis.py` | From the five mask maps: the deterministic dipole-residual bias fields, the two-stage budget test (η_bal ε_res), the product test. | mask maps → `stage3_bias_results.pkl` |
| `sim/stage4_gls.py` | The paper's field estimator: cubic-B-spline model of v_x fitted by generalized least squares to the per-x_c moment samples (heteroscedastic, known σ_k, mild ridge); builds the design matrices and runs the noiseless truncation-bias gate. | maps → `stage4_design.pkl` (190 MB) |
| `sim/stage4_run.py` | Production: GLS field maps and trajectories over the atom-number sweep, matched-filter detection SNR per plane at 10⁴ (0.42 at the first plane, 0.47 pooled), five-sigma atom thresholds, and the dipole bias propagated through the GLS. | design, maps → `stage4_results.pkl` |
| `sim/stage5_alloc.py Na nreal` | Atom allocation across planes: per-plane trajectory-error coefficients by single-plane noise injection, additivity/1/N gates, the optimal allocation and the seed-paired uniform-vs-optimal comparison (m16). | design → `stage5_inject_Na*.pkl` |
| `sim/stage6_fieldmap.py` | The early-plane field-map campaign (planes 0–2): relative/absolute field error vs N_a, total and oscillatory (channelling) matched-filter SNR, the campaign tiers. | design → `stage6_results.pkl` |
| `sim/make_fig4_data.py` | 20 GLS realizations per atom number 10⁴…10⁷: RMS(z) curves, envelope ratios, one representative 10⁷ realization, the noiseless floor. **Regenerating it replaces the anchor.** | design, true_traj → `fig4_data.pkl` |
| `collect_v8_summary.py` | Gathers every headline number into `v8_summary.json` (the master reads this, not the pickles, for most checks). | pickles → `v8_summary.json` |
| `make_fig2.py`, `make_fig3.py` | Closed-form figures (dispersive scaling with detuning; the (F, N_a) regime map). No data files. | — → PDF + check PNG |
| `make_fig4_new.py`, `make_fig5_new.py` | Fig. 4 (needs `stage4_design.pkl`) and Fig. 5 (shipped files suffice). | pickles → PDF + check PNG |
| `make_fig1.py` | Parked script-generated alternative to the schematic; not used by the paper. | — |
| `verify_kinematics.py` | Every kinematic number of Sec. 2 / Table 1 from CODATA constants and the geometry. | — |
| `verify_p2_couplings.py` | The coupling chain (Sec. 3, Tables 3–4) under the 1/e² convention, per-plane budgets, the review-finding closures; run from the root. | `simcode/sim/stage4_results.pkl` |
| `verify_p3_systematics.py` | Every number of Sec. 8 (Coriolis, magnetic, gas, detection pairing, optics budget). | — |
| `verify_g0.py`, `verify_p4_numbers.py` | The correction campaign's decision points; the Secs. 10–11 number set from the result files (run from `simcode_v8/`). | pickles, summary |
| `verify_appB.py` | Appendix B, symbolically (sympy: linearization coefficient, meter commutator, Gaussian-meter pointer shift, kernel identity, Madelung step, moment expansions, narrow-curtain limit) and numerically (free-flight identities and the pointer shift by exact propagation of the joint atom–meter state; the exact decomposition of the verification-grid error; the Sec.-9 partial-trace algebra). | `sim_core`, `v8_config` |
| `verify_glossary.py`, `verify_p8_8.py`, `verify_p8_8b.py` | Appendix C vs body; retired phrases / invariants / bibliography; the final-pass numbers (ramp comparison, single-lobe statistics, N_a exponents, Fig-2 crossover, σ_G, recoil bound, atom numbers). | `work/main.tex`, anchor, `true_traj.pkl` |
| `verify_p7_master.py` | The master: sections A–I recompute and pin every quoted number and run the structural sweeps of `work/main.tex`; section D runs the scripts above as subprocesses and counts their PASS/FAIL lines. Exit 0 only if everything passes. | all of the above |

## 5. Reading the output

- Every verifier prints one line per check: `PASS <name> calc=<recomputed> quoted=<manuscript>` (numbers) or
  `PASS sweep: <name> matches=<n>` (text checks). The last line is the tally. The name tells you where the
  number lives (`Sec 1`, `Table 3`, `Fig-4`, `App B.5`, …).
- A `FAIL` means recomputed ≠ manuscript at the quoted precision (default 2 %). Before concluding anything:
  is the tolerance appropriate for the quoted precision? is a chaotic-regime cell involved (Sec. 2 above)?
  did a regenerated pickle replace the anchor? Only then is it a manuscript error — and then it is real.
- The chain's gates print their own verdicts (`STAGE-1 GATES: 1a PASS 1b PASS`, `GATE 1 worst plane …`,
  `PASS - saved …`); `run_all_v8.sh` stops at the first failure (`set -e`).
- `v8_summary.json` is the human-readable state of the chain: `stage4.snr_1e4` (per-plane SNR at 10⁴),
  `stage4.Na_5sigma` (the 1.43×10⁶ / 8.15×10⁶ / 2.80×10⁷ thresholds), `stage5` (m16), `stage6`, `campaign`
  (tier durations from the flux window), `fig4` (the Fig-4/5 statistics: the anchor's when it is in place).

## 6. Troubleshooting

| symptom | cause / fix |
|---|---|
| `FATAL: stale sim/*.pkl present` | the chain's preflight; `rm simcode_v8/sim/*.pkl` (the shipped small files are regenerated) or use the wrapper, which does this |
| `IndexError` on `sys.argv` | a script that needs arguments was run without them (Spyder F5); see Sec. 3 |
| `ModuleNotFoundError: v8_config` / `No such file: sim/...` | wrong working directory; chain and figure scripts run from `simcode_v8/` |
| `FileNotFoundError: simcode/sim/stage4_results.pkl` | `verify_p2_couplings.py` run from `simcode_v8/`; run it from the root or via the master |
| master FAILs only in the `fig4` / `m16` / `direct pipeline` lines | the anchor was replaced by a regenerated `fig4_data.pkl`; `cp official_fig4_data.pkl simcode_v8/sim/fig4_data.pkl` and re-run `collect_v8_summary.py` |
| a background job stops after a few minutes | the shell killed it on exit; run in the foreground, `tmux`, or stage by stage |
| figures do not appear in Spyder | the scripts use the file backend; open the saved PDF/PNG |
| `verify_appB.py` takes long | ~10 min is normal (2¹⁵ × 512 joint grids); the first stages print progress |
| RAM | < 3 GB peak; stage 4 holds the 190 MB design matrices plus maps |
