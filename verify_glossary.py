#!/usr/bin/env python3
"""verify_glossary.py -- NJP-120031 P8-7g: Appendix C symbol table must not drift from the body.

Assertion 1: every symbol row of the Appendix-C table corresponds to a symbol that occurs in the manuscript
             BODY (everything outside Appendix C). Mechanics: parse the tabular's first column, take each $...$
             token, normalise whitespace, and require a verbatim occurrence in the body.
Assertion 2: the 12 highest-frequency recurring math symbols in the body -- counted mechanically over the
             documented pattern list PATTERNS below -- appear in the table. Universal constants and non-symbol
             words are whitelisted (WHITELIST, printed).
Prints one PASS/FAIL line per check; exit status 1 on any failure. Wired into verify_p7_master.py section D.
"""
import re, sys, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
T = open('work/main.tex').read()
doc = T[T.find(r'\begin{document}'):T.find(r'\begin{thebibliography}')]
iC = doc.find(r'\section{Notation}'); assert iC > 0, "Appendix C not found"
body, appC = doc[:iC], doc[iC:]
tab = appC[appC.find(r'\label{tab:notation}'):appC.find(r'\end{tabular}')]
norm = lambda s: re.sub(r'\s+', '', s)
body_n = norm(body)
npass = nfail = 0
def check(name, ok, detail=''):
    global npass, nfail
    npass += ok; nfail += (not ok)
    print(f"{'PASS' if ok else 'FAIL'} {name:78s} {detail}")

# ---------------- Assertion 1 ----------------
rows = [l for l in tab.split('\n') if ' & ' in l and l.rstrip().endswith(r'\\') and not l.lstrip().startswith(r'\multicolumn')]
symbols = []
for r in rows:
    first = r.split(' & ')[0]
    symbols += re.findall(r'\$([^$]+)\$', first)
missing = [s for s in symbols if norm(s) not in body_n]
check(f"A1 every table symbol occurs verbatim in the body outside App C ({len(symbols)} symbols, {len(rows)} rows)", len(missing) == 0, f"missing: {missing}" if missing else "")

# ---------------- Assertion 2 ----------------
PATTERNS = {  # documented symbol-pattern list (LaTeX-source regexes), shared with the P8-7g inclusion scan
 'x_c': r'x_c\b', 'x_f': r'x_f\b', 'x_gamma': r'x_\\gamma', 'x_a': r'x_a\b', 'z_j': r'z_j\b',
 'w': r'(?<![A-Za-z_\\])w\b(?!_)', 'w_eff': r'w_\\text\{eff\}', 'sigma_0': r'\\sigma_0', 'L': r'(?<![A-Za-z\\])L\b(?!_)', 'T': r'(?<![A-Za-z\\])T\b(?!_)',
 't(z)': r't\(z', 'tau_a': r'\\tau_a', 'v_z': r'v_z', 'v_0': r'v_0\b', 'varphi_1': r'\\varphi_1', 'varphi_gamma': r'\\varphi_\\gamma', 'g_a': r'g_a\b',
 'g_slice': r'g_\\text\{slice\}', 'N_gamma': r'N_\\gamma', 'F': r'(?<![A-Za-z\\])F\b(?!_)', 'Delta': r'\\Delta\b(?!\s*[vxzPp])', 'Gamma': r'\\Gamma\b', 'U_0': r'U_0',
 'tau_1': r'\\tau_1', 'lambda_L': r'\\lambda_L', 'Lambda_eff': r'\\Lambda_\\text\{eff\}', 'Lambda_1': r'\\Lambda_1', 'V': r'(?<![A-Za-z\\])V\b(?!_)',
 'D': r'(?<![A-Za-z\\])D\b(?!_)', 'D_1': r'\\mathcal\{D\}_1', 'P_sc': r'P_\\text\{sc\}', 'N_a': r'N_a\b', 'sigma_p': r'\\sigma_p', 'sigma_x': r'\\sigma_x',
 'G': r'(?<![A-Za-z\\])G\b(?!_)', 'K': r'(?<![A-Za-z\\])K\b(?!_)', 'rho_w': r'\\rho_w', 'J_w': r'J_w', 'scriptA': r'\\mathcal\{A\}', 'v_x': r'v_x\b',
 'ftilde': r'\\tilde\{f\}', 'Ahat': r'\\hat\{A\}', 'weakval': r'\\weakval', 'psi': r'\\psi\b', 'rho': r'\\rho\b(?!_)', 'J': r'(?<![A-Za-z\\_])J\b(?!_)',
 'delta_p': r'\\delta p', 'Phi_0': r'\\Phi_0', 'eta_bal': r'\\eta_\\text\{bal\}', 'eps_res': r'\\varepsilon_\\text\{res\}', 'V_dipole': r'V_\\text\{dipole\}',
 'p_gamma': r'p_\\gamma', 'hbar': r'\\hbar', 'M': r'(?<![A-Za-z\\])M\b(?!_)', 'g': r'(?<![A-Za-z\\_])g\b(?!_)', 'SNR': r'\bSNR\b',
}
WHITELIST = {'hbar': 'universal constant', 'M': 'atomic mass (universal constant of the problem)', 'g': 'gravitational acceleration', 'SNR': 'a word, not a symbol'}
TABLE_TOKEN = {'x_c': 'x_c', 'x_f': 'x_f', 'x_gamma': r'x_\gamma', 'x_a': 'x_a', 'z_j': 'z_j', 'w': '$w$', 'w_eff': r'w_\text{eff}', 'sigma_0': r'\sigma_0', 'L': '$L$', 'T': '$T$',
 't(z)': 't(z)', 'tau_a': r'\tau_a', 'v_z': 'v_z(z)', 'v_0': 'v_0', 'varphi_1': r'\varphi_1', 'varphi_gamma': r'\varphi_\gamma', 'g_a': 'g_a', 'g_slice': r'g_\text{slice}',
 'N_gamma': r'N_\gamma^{(1)}', 'F': '$F$', 'Delta': r'\Delta', 'Gamma': r'\Gamma', 'U_0': 'U_0', 'tau_1': r'\tau_1', 'lambda_L': r'\lambda_L', 'Lambda_eff': r'\Lambda_\text{eff}',
 'Lambda_1': r'\Lambda_1', 'V': '$V$', 'D': '$D$', 'D_1': r'\mathcal{D}_1', 'P_sc': r'P_\text{sc}', 'N_a': 'N_a', 'sigma_p': r'\sigma_p', 'sigma_x': r'\sigma_x', 'G': 'G(x_c',
 'K': 'K(x_f', 'rho_w': r'\rho_w', 'J_w': 'J_w', 'scriptA': r'\mathcal{A}', 'v_x': 'v_x', 'ftilde': r'\tilde{f}', 'Ahat': r'\hat{A}', 'weakval': r'\weakval', 'psi': r'\psi',
 'rho': r'\rho$', 'J': '$J$', 'delta_p': r'\delta p_\gamma', 'Phi_0': r'\Phi_0', 'eta_bal': r'\eta_\text{bal}', 'eps_res': r'\varepsilon_\text{res}', 'V_dipole': r'V_\text{dipole}', 'p_gamma': r'\hat{p}_\gamma'}
counts = {k: len(re.findall(p, body)) for k, p in PATTERNS.items()}
ranked = [k for k, _ in sorted(counts.items(), key=lambda kv: -kv[1]) if k not in WHITELIST]
top12 = ranked[:12]
print("      whitelist (excluded from the top-12 test): " + ", ".join(f"{k} [{v}]" for k, v in WHITELIST.items()))
print("      top-12 recurring symbols by body occurrence: " + ", ".join(f"{k}({counts[k]})" for k in top12))
absent = [k for k in top12 if TABLE_TOKEN[k] not in tab]
check("A2 the 12 highest-frequency recurring symbols of the body appear in the table", len(absent) == 0, f"absent: {absent}" if absent else "")
# bonus: inclusion-rule coverage -- every non-whitelisted pattern present in >= 2 sections is in the table
secs = re.split(r'\\section\{', body)
cov_missing = [k for k, p in PATTERNS.items() if k not in WHITELIST and sum(bool(re.search(p, s)) for s in secs) >= 2 and TABLE_TOKEN[k] not in tab]
check("A3 inclusion rule: every documented symbol used in >= 2 body sections is in the table", len(cov_missing) == 0, f"missing: {cov_missing}" if cov_missing else "")
check("A4 the Sec-2 notation paragraph points to Appendix C", r'Appendix~\ref{app:notation}' in body)
print(f"GLOSSARY VERIFICATION: {npass} checks passed, {nfail} failed")
sys.exit(0 if nfail == 0 else 1)
