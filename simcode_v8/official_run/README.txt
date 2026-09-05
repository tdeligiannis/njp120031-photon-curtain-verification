official_run/v8_run_log.txt -- console log of the authors' first chain run on the official machine (2026-07-04,
Python 3.12.11 / numpy 2.4.0 / scipy 1.18.0). PROVENANCE NOTES (harness audit, 2026-08-25):
- This log ends at a stage-5 crash (a bug fixed in the shipped stage5_alloc.py); the completed run that followed
  (resume_v8.sh) produced the Fig-4 data set shipped as official_fig4_data.pkl, but its log is not available.
- The stage-2 line of this log (final trajectory RMS 13817.7 +- 11710.4 um) is a DIFFERENT set of 30 noise
  realizations of the same chaotic-regime cell than the shipped sim/stage2_results.pkl (12780.4 +- 11397.6 um,
  n = 30); the two agree within the standard error of the mean (~2.1 mm). The manuscript, Fig. 5, v8_summary.json
  and the master verification all quote the SHIPPED run (~13 mm). Likewise the allocation ratio at 1e6 quoted in
  the manuscript (1.51 +- 0.31) is the shipped run's (the crash preceded that stage in this log).
- Every deterministic quantity in this log (gates, couplings, SNRs, field errors, dipole biases) reproduces
  bit-for-bit from the shipped code on other machines.
