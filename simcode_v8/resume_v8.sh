#!/usr/bin/env bash
# NJP-120031 v8.1 RESUME after the stage-5 m16 fix. Run from simcode_v8/ on the machine that
# already completed the chain through stage4 (maps, design, stage4_results on disk).
# 1) Drop the fixed sim/stage5_alloc.py in place first (run_alloc def moved above __main__).
# Resume is fast: stage5 Na=1e6 skips its finished A-loop and runs only the m16 block.
set -euo pipefail
LOG=v8_run_log.txt
{ echo; echo "== v8.1 RESUME $(date -u +%FT%TZ) (stage5 m16 fix) =="
  T0=$SECONDS
  step(){ echo; echo "---- [$((SECONDS-T0))s] $* ----"; }
  step "stage5_alloc Na=1e6 (m16 only)"; python3 sim/stage5_alloc.py 1e6 12
  step "stage5_alloc Na=1e7";            python3 sim/stage5_alloc.py 1e7 12
  step "stage6_fieldmap";                python3 sim/stage6_fieldmap.py
  step "make_fig4_data";                 python3 sim/make_fig4_data.py
  step "collect summary";                python3 collect_v8_summary.py
  echo; echo "== resume complete in $((SECONDS-T0)) s =="
} 2>&1 | tee -a "$LOG"
