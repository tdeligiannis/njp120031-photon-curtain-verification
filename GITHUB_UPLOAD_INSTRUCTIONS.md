# Publishing the NJP-120031 code: GitHub + Zenodo DOI

You have `NJP120031_code_repo_ready_to_push.zip`: a git repository (one commit, branch `main`) containing the
supplementary code, the verification harness, the official Fig-4 anchor, the small result files, `README.md`,
`LOCAL_RUN_GUIDE.md`, an MIT `LICENSE`, `CITATION.cff`, `.zenodo.json`, `requirements.txt` and a `.gitignore`
that excludes the regenerated large intermediates and figure outputs.

## 0. Timing — read first

The harness checks `work/main.tex`, so the **released** version of the repository must contain the manuscript as
actually submitted. Sequence:

1. Now: create the GitHub repository and push the current commit (private or public — your choice; Zenodo needs it
   public at release time).
2. After Jessica's fixes and the release gate (master 178/0 and letter verifier 7/7 on the final files): replace
   `work/main.tex` with the final manuscript, update `verification_log_master.txt` (and README counts if they change),
   commit, then **tag `v1.0.0` and publish a GitHub release**.
3. DOI: either let the Zenodo–GitHub integration mint it at the release (Sec. 3A), or reserve a DOI on Zenodo first so
   you can write it into the Data availability statement before submitting (Sec. 3B — recommended, because the
   manuscript needs the DOI string and the release should contain the manuscript that contains the DOI).

## 1. Terminal path (Linux/macOS, ~10 minutes)

```bash
unzip NJP120031_code_repo_ready_to_push.zip && cd repo

# your git identity (the initial commit carries a placeholder author; fix it once)
git config user.name  "Theodore Deligiannis"
git config user.email "<your-email>"
git commit --amend --reset-author --no-edit

# sanity: the harness runs from the repo root on shipped files alone (~10 min); skip if short of time
python3 verify_glossary.py && python3 verify_p8_8.py && python3 verify_p8_8b.py
# python3 verify_p7_master.py     # full master, 178 checks

# create the GitHub repository (either of the two)
gh auth status || gh auth login                     # GitHub CLI, one-time
gh repo create <ORG-OR-USER>/njp120031-photon-curtain-verification --public --source=. --remote=origin \
    --description "Supplementary code and verification harness for NJP-120031 (photon-curtain weak measurement)"
# -- or, without gh: create an empty repo in the web UI, then:
# git remote add origin git@github.com:<ORG-OR-USER>/njp120031-photon-curtain-verification.git

git push -u origin main

# fill in the repository URL in CITATION.cff (placeholder <ORG-OR-USER>), commit, push
sed -i 's#<ORG-OR-USER>#<your-org-or-user>#' CITATION.cff && git commit -am "CITATION: repository URL" && git push
```

Release (after step 0.2):
```bash
git tag -a v1.0.0 -m "Code verified against the resubmitted manuscript (NJP-120031, September 2026)"
git push origin v1.0.0
gh release create v1.0.0 --title "v1.0.0 — resubmission" \
   --notes "Verification harness: 178 checks, all live, all passing against the resubmitted manuscript (see verification_log_master.txt)."
```

## 2. What to check in the GitHub UI after the push
- The README renders (the package README is the front page); `LOCAL_RUN_GUIDE.md` is linked from it.
- `official_fig4_data.pkl` (176 KB) and the ten small `simcode_v8/sim/*.pkl` are present; no `stage2_maps.pkl`,
  `stage3_maps_*.pkl` or `stage4_design.pkl` (excluded by `.gitignore`).
- Settings → General → "Include Git LFS objects in archives" is irrelevant (no LFS); nothing else to configure.
- Optional: enable "Discussions" off, "Issues" on (referees may file questions there).

## 3. Minting the DOI

### 3A. Zenodo–GitHub integration (DOI created at the release)
1. Log in to zenodo.org with GitHub → "GitHub" in the account menu → flip the switch for the repository.
2. Publish the GitHub release `v1.0.0` (Sec. 1). Zenodo archives the release automatically and mints a DOI; the
   `.zenodo.json` in the repo supplies title, authors, license and description.
3. Copy the **concept DOI** (resolves to the latest version) into the Data availability statement; the version DOI
   identifies exactly the released files.

### 3B. Reserve the DOI first (recommended for this submission)
1. zenodo.org → New upload → "Reserve DOI" (button next to the DOI field). You get the DOI string immediately, before
   publishing anything.
2. Write that DOI into `main.tex` (Data availability: `[DOI to be inserted before publication]` → the DOI), rebuild,
   re-run the harness, commit to the repository, tag and release `v1.0.0`.
3. On the reserved Zenodo record: upload the release zip from GitHub (or use the integration, which will link to the
   same concept), fill authors/affiliations/license (MIT), related identifier "is supplement to" the article, and
   **publish**. The DOI then resolves to the exact code the manuscript cites.

## 4. Prompt for a fresh Claude Code (or Cowork) session

Paste this in a session opened in the directory where you unzipped the repo, with `gh` authenticated:

> You are in a git repository (branch `main`, one commit) containing the supplementary code and verification
> harness for our NJP manuscript NJP-120031. Read `README.md` and `GITHUB_UPLOAD_INSTRUCTIONS.md` first.
> Tasks, in order, stopping to show me the result of each before continuing:
> 1. Set the git identity to Theodore Deligiannis / <email> and amend the initial commit's author
>    (`git commit --amend --reset-author --no-edit`). Do not change any file content.
> 2. Run the three quick verifiers (`verify_glossary.py`, `verify_p8_8.py`, `verify_p8_8b.py`) in the foreground and
>    report their last lines. Do not edit any verifier or any file under `simcode_v8/` for any reason; if something
>    fails, stop and report.
> 3. Create the public GitHub repository `<ORG-OR-USER>/njp120031-photon-curtain-verification` with `gh repo create`
>    using this directory as source and `origin` as remote, push `main`, and give me the URL.
> 4. Replace the placeholder `<ORG-OR-USER>` in `CITATION.cff` with the actual owner, commit ("CITATION: repository
>    URL") and push.
> 5. Do NOT create a tag or a release, and do NOT touch Zenodo — the release happens only after the final manuscript
>    copy is committed; I will tell you when.
> 6. Finally, list every file that was pushed (`git ls-files`) and confirm that no `stage2_maps.pkl`,
>    `stage3_maps_*.pkl` or `stage4_design.pkl` is tracked.
> Never run anything in the background; never force-push; ask before any action not listed above.

For the later release step, a second prompt:

> The final manuscript copy is now in `work/main.tex` and `verification_log_master.txt` is updated. Confirm
> `python3 verify_p7_master.py` ends with `178 PASS / 0 FAIL` (foreground, ~10 min). Then commit ("Final manuscript
> copy for release"), tag `v1.0.0` with the message "Code verified against the resubmitted manuscript (NJP-120031,
> September 2026)", push the tag, and create the GitHub release with `gh release create` using the notes in
> `GITHUB_UPLOAD_INSTRUCTIONS.md` Sec. 1. Report the release URL. Do nothing else.

## 5. What NOT to put in the repository
- The referee reports, the decision letter, the response letter, STATE.md, the edit scripts and the campaign records
  (internal); the response and highlighted PDFs go to the journal, not to the code archive.
- API keys, credentials, or any file from your home directory.
