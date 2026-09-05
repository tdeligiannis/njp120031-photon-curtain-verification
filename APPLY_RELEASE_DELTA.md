# Release delta for tdeligiannis/njp120031-photon-curtain-verification (apply BEFORE tagging v1.0.0)

Files in this folder replace their namesakes at the repository root:
  README.md           new title ("following Kocsis–Steinberg"), September 2026, 179-check statements, do-not-delete note on work/main.tex
  CITATION.cff        full paper title (Zenodo inherits DOI metadata from CITATION/.zenodo), repository URL
  .zenodo.json        full paper title; 179 checks; related identifier
  verify_p7_master.py IOP lowercase-style sweeps (figure/appendix/Reading figure); A1 sweep by intent
  verify_glossary.py  case-insensitive Appendix-C pointer check
  verify_p8_8.py      35 references; m16 string tolerated twice (Sec 10.2 + App D.5); D.5 seed-pair pin (12, never 20)

Then, and only then, the manuscript copy:
  work/main.tex       <- the FINAL main.tex (the one that compiles to the 38-page PDF of 5 Sept, with '12 seed pairs',
                         eq:weak_value_def, British spelling, the ", The relevant coherence" fix and the abstract-count
                         wording settled). NEVER delete this file: three verifiers sweep it.

Sequence:
  cp <delta>/* .   && cp <delta>/.zenodo.json .
  cp <final>/main.tex work/main.tex
  python3 verify_glossary.py && python3 verify_p8_8.py && python3 verify_p8_8b.py     # quick (expect 4/0, 56/0, 24/0)
  python3 verify_p7_master.py | tee verification_log_master.txt                        # ~10 min; expect 179 PASS / 0 FAIL
  python3 - <<'PY'                                                                     # AST scan: 0 literal-value checks
import ast, glob
files=['verify_p7_master.py','verify_appB.py','verify_glossary.py','verify_p8_8.py','verify_p8_8b.py']+sorted(glob.glob('simcode_v8/verify_*.py'))
h=0
for f in files:
    src=open(f).read(); tree=ast.parse(src); helpers=set()
    for n in ast.walk(tree):
        if isinstance(n, ast.FunctionDef) and len(n.args.args)>=2:
            a=[x.arg for x in n.args.args[:2]]
            if a[0] in ('name','label','tag') and a[1] in ('val','value','calc','x','v'): helpers.add(n.name)
    for n in ast.walk(tree):
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id in helpers and len(n.args)>=2 and isinstance(n.args[1], ast.Constant) and isinstance(n.args[1].value,(int,float)): h+=1
print("literal-value checks:", h)
PY
  git add -A && git commit -m "Release prep: final manuscript copy, harness for the resubmission layout, title/date/metadata" && git push
  git tag -a v1.0.0 -m "Code verified against the resubmitted manuscript (NJP-120031, September 2026)" && git push origin v1.0.0
  gh release create v1.0.0 --title "v1.0.0 — resubmission" --notes "Verification harness: 179 checks, all live, all passing against the resubmitted manuscript (verification_log_master.txt)."

Zenodo: reserve the DOI first, write it into main.tex's Data availability statement, rebuild, re-run the master
(the copy in work/ must be the version WITH the DOI), commit, then tag/release/publish. Keep the code record separate
from the April preprint record; link them with related identifiers (code: isSupplementTo the article/preprint DOI;
preprint v2: isSupplementedBy the code DOI).
