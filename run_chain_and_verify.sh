#!/usr/bin/env bash
# NJP-120031 supplementary code -- full reproduction with the official Fig-4 anchor preserved.
# Run from supplementary_code/ (the directory containing this script). Needs python3 + numpy/scipy/sympy/mpmath/matplotlib.
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"
echo "== NJP-120031 full reproduction $(date -u +%FT%TZ) =="
echo "manuscript copy md5: $(md5sum work/main.tex | cut -c1-32)   (compare with the Overleaf root main.tex)"
[ -f official_fig4_data.pkl ] || { echo "FATAL: official_fig4_data.pkl missing"; exit 1; }
cd simcode_v8
rm -f sim/*.pkl                                   # the chain's preflight refuses stale pickles; the small shipped
                                                  # result files are regenerated (deterministic ones identically)
bash run_all_v8.sh                                # ~8 min (workstation) to ~40 min (laptop); ~700 MB of intermediates
cp ../official_fig4_data.pkl sim/fig4_data.pkl    # re-install the official anchor: make_fig4_data.py regenerated a
                                                  # platform-variant set (chaotic regime); the manuscript quotes the official one
python3 collect_v8_summary.py                     # summary with the anchor in place
cd ..
echo "== master verification =="
python3 verify_p7_master.py
