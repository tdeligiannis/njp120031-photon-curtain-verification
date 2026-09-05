<!-- Repository header (GitHub); the package README follows unchanged. -->
**Repository:** supplementary code and verification harness for NJP-120031 (Fabianiak & Deligiannis, New Journal of Physics).
**Status:** version 1.0.0 = the code verified against the resubmitted manuscript. **License:** MIT (LICENSE). **Cite:** CITATION.cff; archived DOI via Zenodo (see the manuscript's Data availability statement).
**Quick start:** `pip install -r requirements.txt`, then Sec. 3 below; the full reproduction is `bash run_chain_and_verify.sh` (Sec. 4).
**Note on `work/main.tex`:** a verification copy of the manuscript, included because the harness checks the text (do not delete it: `verify_p7_master.py`, `verify_glossary.py` and `verify_p8_8.py` sweep it); at the tagged release it is byte-identical to the journal submission.

# NJP-120031 — Supplementary code and data (revised manuscript, September 2026)

Fabianiak & Deligiannis, "Weak-measurement reconstruction of the conditional momentum field for atomic
matter waves: a dispersive photon-curtain protocol following Kocsis–Steinberg", New Journal of Physics.

This folder contains the complete Python source that generates Figures 2–5 and **every simulation- or
formula-derived number quoted in the manuscript**, together with the verification harness that re-derives
each of those numbers and checks it against the manuscript text, and the small result files needed to
reproduce the figures without re-running the simulation chain.

Everything is deterministic: all random number generators are explicitly seeded. Three quantities live in a
chaotic regime of the estimator (the direct-pipeline trajectory RMS at 10^4 atoms per point, the
allocation-gain ratio at 10^6, and the Fig-4 trajectory statistics) and vary at the few-percent level
between machines. The manuscript quotes, for each of them, the value of a SHIPPED artifact: the Fig-4
statistics come from the authors' official run (shipped as the anchor `official_fig4_data.pkl`), the other
two from the shipped result files `sim/stage2_results.pkl` and `sim/stage5_inject_Na1e6.pkl` (see "Provenance").

## 1. Layout

    supplementary_code/
    ├── README.md                      this file
    ├── verify_p7_master.py            MASTER VERIFICATION (179 checks): every quoted number and every
    │                                  structural invariant of main.tex; runs the scripts below as subprocesses
    ├── verify_appB.py                 Appendix B: every identity checked symbolically (sympy) and numerically
    ├── verify_glossary.py             Appendix C symbol table vs the manuscript body
    ├── verify_p8_8.py                 stale-pattern sweep, do-not-regress invariants, bibliography cross-checks
    ├── verify_p8_8b.py                pins for the numbers introduced in the final revision pass
    ├── verification_log_master.txt    the master's output for the shipped manuscript (179 PASS / 0 FAIL)
    ├── run_chain_and_verify.sh        full reproduction: simulation chain -> anchor re-installed -> master
    ├── official_fig4_data.pkl         the official Fig-4 data set (anchor; see Provenance)
    ├── work/main.tex                  verification copy of the manuscript (identical to ../main.tex at packaging)
    ├── simcode/sim/stage4_results.pkl v7-era result file read by ONE closure check in verify_p2_couplings.py
    └── simcode_v8/                    the simulation chain
        ├── sim/                       chain modules (v8_config, sim_core, stage2_* ... stage6_*, make_fig4_data)
        │   ├── fig4_data.pkl          = official_fig4_data.pkl (anchor copy, read by the figure scripts)
        │   └── *.pkl (small)          shipped result files: stage1_verification, stage2ref_gate, stage2_results,
        │                              stage3_bias_results, stage4_results, stage5_inject_Na1e6/1e7,
        │                              stage6_results, true_traj  (the large intermediates are NOT shipped)
        ├── run_all_v8.sh              the official chain script (as run for the record; regenerates fig4_data.pkl)
        ├── resume_v8.sh               historical resume script (provenance only)
        ├── collect_v8_summary.py      collects every headline number into v8_summary.json
        ├── v8_summary.json            summary of the shipped chain state (fig4 section = the anchor)
        ├── make_fig2.py ... make_fig5_new.py   figure generators (Figs. 2-5); make_fig1.py is a parked
        │                              script-generated alternative to the schematic (not used)
        ├── verify_g0.py, verify_kinematics.py, verify_p2_couplings.py, verify_p3_systematics.py,
        │   verify_p4_numbers.py       gate verifiers of the correction campaign (subprocesses of the master)
        ├── official_run/              the authors' official run log (+ README on its provenance)
        ├── README_v8.md               chain notes from the correction campaign (historical)
        └── legacy/                    superseded figure scripts kept for provenance (see Sec. 6)

## 2. Environment

Python >= 3.10 with numpy, scipy, matplotlib, sympy and mpmath. Reference environment of the shipped
results: Python 3.12.3, numpy 2.4.4, scipy 1.17.1, matplotlib 3.10.8, sympy 1.14, mpmath 1.3. No other
dependencies; no compiled extensions. Run everything from the directory indicated below.

## 3. Quick start (< 2 minutes; uses the shipped result files)

    cd supplementary_code
    python3 verify_glossary.py          # 4 checks    (~1 s)
    python3 verify_p8_8.py              # 51 checks   (~1 s)
    python3 verify_p8_8b.py             # 24 checks   (~1 min)
    python3 verify_p7_master.py         # the full master, 179 checks (~10 min; the shipped files suffice)
    cd simcode_v8
    python3 make_fig2.py                # Fig. 2  (closed forms; no data files)
    python3 make_fig3.py                # Fig. 3  (closed forms; no data files)
    python3 make_fig5_new.py            # Fig. 5  (reads sim/fig4_data.pkl, stage6_results.pkl, stage4_results.pkl)
    python3 verify_p4_numbers.py        # numbers of Secs. 10-11 from the shipped result files

Figure 4 additionally needs sim/stage4_design.pkl (190 MB, not shipped): regenerate it with the chain (Sec. 4)
and then run `python3 make_fig4_new.py`. Regenerated figures reproduce the pinned numbers on any stack
and are pixel-identical to the shipped PDFs under the reference matplotlib (3.10.8); newer matplotlib
versions change bytes, not content.

## 4. Full reproduction (simulation chain + master verification; ~10-40 min)

    cd supplementary_code
    bash run_chain_and_verify.sh

This (i) clears sim/*.pkl (the chain's preflight refuses to start over stale pickles), (ii) runs the
official chain script `simcode_v8/run_all_v8.sh` — 20 curtain-plane record maps, 30 direct-pipeline noise
realizations, five dipole-residual mask maps, the GLS design and runs, the allocation study, the field
maps, the Fig-4 trajectory data and the summary; ~8 min on a workstation, 20-40 min on a laptop; ~700 MB
of intermediates — (iii) re-installs the official Fig-4 anchor (see Provenance) and re-collects the
summary, and (iv) runs the master verification, which must report `179 PASS / 0 FAIL`.

Gates inside the chain abort loudly if violated: v8_config self-checks; the Sec.-5.4 24-point grid
(mean 1.04 %, max 2.52 %); the free-fall gravity-map identity (max 1.64 %); the cumulative-vs-direct
estimator identity (worst plane 0.16 %). All of these are deterministic and reproduce exactly.

If you run `run_all_v8.sh` on its own, note that its `make_fig4_data` step regenerates `sim/fig4_data.pkl`;
the chaotic-regime cells then differ from the manuscript at the few-percent level (the shipped
`v8_summary.json` records the official values), which is expected. `run_chain_and_verify.sh` restores the
anchor for exactly this reason.

Worked example (packaging test, 2026-08-25, on a machine different from the official run): the chain rebuilt
from this folder reproduced every deterministic gate exactly (Stage 1a 1.04/2.52 %; 1b 1.64 %; identity
0.16 %; GLS gate; SNR 0.42/0.47; m16 1.113 at 10^7), the regenerated Fig-4 data gave 1237 um at 10^4 and
183.9 um at 10^7 against the anchor's 1252.8 and 183.8 um (the chaotic-regime spread), and the master
verification returned 179 PASS / 0 FAIL after the anchor was re-installed.

Note on long runs: if your shell kills background jobs on exit, run the wrapper in a foreground terminal or
run the stages of `run_all_v8.sh` one by one; every stage checkpoints its own pickle.

## 5. Verification harness

`verify_p7_master.py` recomputes every number quoted in the manuscript from closed forms, the production
result files or the official run, and PASS/FAILs each against the quoted value at the quoted precision;
its later sections run structural sweeps of `work/main.tex` (retired phrases must be absent, invariants
present, every symbol of the Appendix-C table used in the body, etc.) and call the subprocess verifiers.
The output for the shipped manuscript is `verification_log_master.txt`. Each subprocess verifier can be run
on its own from this directory (they locate their inputs relative to their own location):

| script | what it establishes |
|---|---|
| verify_kinematics.py | Sec. 2 and Table 1: free-fall kinematics, per-plane times, fringe spacing (26 checks) |
| verify_p2_couplings.py | Sec. 3 and Tables 3-4: the coupling chain under the 1/e^2 convention, per-plane budgets (69) |
| verify_p3_systematics.py | Sec. 8: Coriolis, magnetic, background-gas, detection and optics budgets (48) |
| verify_g0.py, verify_p4_numbers.py | independent rederivation of the design decisions; the Secs. 10-11 number set |
| verify_appB.py | Appendix B: linearization coefficient, meter commutator, Gaussian-meter pointer shift, kernel identity, Madelung step, moment expansions, narrow-curtain limit (sympy); free-flight identities and the pointer shift on the exact two-Gaussian state by joint atom-meter propagation; the exact decomposition of the verification-grid error; the Sec.-9 partial-trace algebra (35 checks, ~10 min) |
| verify_glossary.py | Appendix C: every table symbol occurs in the body; the top-12 body symbols are tabulated (4) |
| verify_p8_8.py | 22 retired phrases absent; the do-not-regress invariants present; bibliography and label cross-checks (51) |
| verify_p8_8b.py | the numbers of the final revision pass: far-field ramp comparison, single-lobe pattern statistics, RMS-vs-N_a exponents, Fig-2 crossover, diffracted-slice width, recoil bound, atom numbers (24) |

## 6. Provenance and notes

- **Official anchor.** `official_fig4_data.pkl` (and its copy `simcode_v8/sim/fig4_data.pkl`) is the Fig-4
  trajectory data of the authors' official chain run (2026-07-04): 20 noise realizations at each of
  10^4 ... 10^7 atoms per point plus the noiseless-pipeline floor. The manuscript's Fig-4/Fig-5 trajectory
  statistics (1.25 +- 0.84 mm at 10^4; 184 +- 93 um at 10^7; floor 3.4 um) are pinned to it. The official
  machine's (partial) run log is in `simcode_v8/official_run/` with a README on its provenance; the two other
  chaotic-regime quantities (direct-pipeline RMS at 1e4, allocation ratio at 1e6) are quoted from the shipped
  result files, which are the run of record for them, and the master checks them live.
- **Naming.** The manuscript calls the recommended apparatus "the operating configuration"; in the code it is
  still labelled `Configuration E` (constants block `v8_config.py`, docstrings). The two are the same thing.
  The single-pass comparison case of Fig. 2 and Table 3 ("single-pass free-space reference", 780 nm,
  462 uW) is the code's `Config-D` reference block.
- **Figure 2 (b).** The dipole-kick curve is drawn with the 1/e intensity radius w_eff, as in Eq. (20) of the
  manuscript; an earlier version of `make_fig2.py` used the 1/e^2 radius w (curve too low by sqrt 2,
  crossover 0.50 GHz instead of 0.35 GHz). The superseded script is kept in `legacy/P6_fig2_w-slip/`.
- **Legacy folder.** `legacy/stale_v7_fig_scripts/` holds the figure scripts of the pre-correction (v7)
  bundle; `legacy/build_simcode_v8.py` is the builder that produced this chain from the v7 code, and
  `legacy/NJP120031_v7_final_verification_output.txt` is the v7 verification record. None of these is needed
  to reproduce the manuscript.
- **v7 input.** `simcode/sim/stage4_results.pkl` is a result file of the pre-correction chain, read by one
  block of `verify_p2_couplings.py` that closes a review finding against the earlier numbers. It is not used
  by any figure or by any current manuscript number.
- **Manuscript copy.** `work/main.tex` is the manuscript as verified; it is byte-identical to the Overleaf
  root `main.tex` of this bundle at packaging time (the md5 is printed by `run_chain_and_verify.sh`).
