#!/usr/bin/env bash
# NJP-120031 v8 full simulation chain. Run from simcode_v8/. ~20-40 min on a modern laptop.
# Produces sim/*.pkl + v8_run_log.txt + v8_summary.json. Send back: v8_summary.json,
# v8_run_log.txt, and the sim/ pickles EXCEPT stage2_maps.pkl / stage3_maps_*.pkl /
# stage4_design.pkl (large intermediates; keep locally in case of follow-ups).
set -euo pipefail
LOG=v8_run_log.txt
{ echo "== NJP-120031 v8 chain start $(date -u +%FT%TZ) =="
  python3 --version; python3 -c "import numpy,scipy;print('numpy',numpy.__version__,'scipy',scipy.__version__)"
  echo "== preflight =="
  echo "chain version: v8.1 (six-fix forensic build, 2026-07-04)"
  echo "running from: $(readlink -f .)"
  case "$(readlink -f .)" in *Trash*|*trash*)
    echo "FATAL: running from the Trash. Your shell followed a deleted directory."
    echo "       Open a fresh terminal, cd to the real extraction, and rerun."; exit 1;; esac
  grep -q GRID_REL sim/stage2_maps.py    || { echo "FATAL: v8.0 tree detected (f2 fix absent). Extract the v8.1 zip."; exit 1; }
  grep -q "side='right'" sim/stage2_maps.py || { echo "FATAL: f3 fix absent. Extract the v8.1 zip."; exit 1; }
  grep -q cum_integrate sim/sim_core.py  || { echo "FATAL: f4 fix absent. Extract the v8.1 zip."; exit 1; }
  grep -q a_w_inc sim/stage2_maps.py     || { echo "FATAL: f5 fix absent. Extract the v8.1 zip."; exit 1; }
  if ls sim/*.pkl >/dev/null 2>&1; then
    echo "FATAL: stale sim/*.pkl present (resume logic would silently reuse them)."
    echo "       Run: rm sim/*.pkl   -- or extract a fresh tree."; exit 1; fi
  echo "preflight OK: v8.1 tree, clean pickle state"
  T0=$SECONDS
  step(){ echo; echo "---- [$((SECONDS-T0))s] $* ----"; }
  step "v8_config self-checks";        python3 sim/v8_config.py
  step "stage 1 gates (sim_core)";     python3 sim/sim_core.py
  step "stage2_maps (20 planes)";      python3 sim/stage2_maps.py 0 20
  step "stage2_reference gate";        python3 sim/stage2_reference.py
  step "stage2_run (30 realizations)"; python3 sim/stage2_run.py 30
  for MI in 0 1 2 3 4; do
    step "stage3_maps mask $MI";       python3 sim/stage3_maps.py $MI 0 20
  done
  step "stage3_analysis";              python3 sim/stage3_analysis.py
  step "stage4_gls (design + gate)";   python3 sim/stage4_gls.py
  step "stage4_run";                   python3 sim/stage4_run.py
  step "stage5_alloc Na=1e6";          python3 sim/stage5_alloc.py 1e6 12
  step "stage5_alloc Na=1e7";          python3 sim/stage5_alloc.py 1e7 12
  step "stage6_fieldmap";              python3 sim/stage6_fieldmap.py
  step "make_fig4_data";               python3 sim/make_fig4_data.py
  step "collect summary";              python3 collect_v8_summary.py
  echo; echo "== chain complete in $((SECONDS-T0)) s =="
} 2>&1 | tee "$LOG"
