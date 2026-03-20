# JOSS Submission Checklist

Pre-submission checklist for the Journal of Open Source Software (JOSS).
JOSS review criteria: https://joss.readthedocs.io/en/latest/review_criteria.html
JOSS paper format: https://joss.readthedocs.io/en/latest/paper.html (750-1,750 words)

---

## Section A: JOSS Requirements Compliance

### Software Requirements

| # | Requirement | Status | Evidence |
|---|-------------|--------|----------|
| 1 | Open source license (OSI-approved) | MET | BSD-3-Clause in `LICENSE`; declared in `CITATION.cff` |
| 2 | Version control (public repository) | MET | https://github.com/USACE-RMC/Numerics |
| 3 | README with project description | MET | `README.md` with description, install, docs, badges |
| 4 | Installation instructions | MET | NuGet install in `README.md` and `docs/getting-started.md` |
| 5 | Usage examples | MET | `docs/getting-started.md`, `docs/index.md` Quick Start, paper Example section |
| 6 | API documentation | MET | 25+ Markdown files in `docs/`; XML documentation generated from source |
| 7 | Community guidelines (CONTRIBUTING) | MET | `CONTRIBUTING.md` with bug reports, PRs, DCO, security policy |
| 8 | Code of Conduct | MET | `CODE_OF_CONDUCT.md` (Contributor Covenant v2.1) |
| 9 | Automated tests | MET | 1,006 MSTest methods across 150 test classes in `Test_Numerics/` |
| 10 | Continuous integration | MET | 3 GitHub Actions workflows (`Integration.yml`, `Snapshot.yml`, `Release.yml`) |
| 11 | Substantial scholarly effort | MET | 60,000+ LOC library, 34,000+ LOC tests, 248 source files |
| 12 | Sustained development history | MET | 2.5 years (Sep 2023 - Mar 2026), 270+ commits, 7+ contributors |
| 13 | Statement of need | MET | Paper includes Statement of Need section |
| 14 | State of the field / related work | MET | Paper includes State of the Field section |
| 15 | Zero or documented dependencies | MET | Zero runtime dependencies for .NET 8+; polyfills only for .NET Framework 4.8.1 |
| 16 | Functionality matches claims | MET | All claims verified: 43 distributions, 8 MCMC samplers, 5+ optimizers, copulas, ML, bootstrap |
| 17 | Software installable | REMAINING | Publish `RMC.Numerics` to nuget.org (see Step 4 below) |

### Paper Requirements

| # | Requirement | Status | Evidence |
|---|-------------|--------|----------|
| 18 | Paper in `paper/paper.md` | MET | `paper/paper.md` with YAML frontmatter |
| 19 | Bibliography in `paper/paper.bib` | MET | `paper/paper.bib` with 29 references |
| 20 | Summary section | MET | Lines 39-41 |
| 21 | Statement of Need section | MET | Lines 43-54 |
| 22 | State of the Field section | MET | Lines 56-62 |
| 23 | Research impact / ongoing projects | MET | Lines 64-75 (6 production applications) |
| 24 | Acknowledgements | MET | Lines 119-121 |
| 25 | AI Usage Disclosure | MET | Lines 115-117 |
| 26 | Word count 750-1,750 | MET | ~1,215 words (within range) |
| 27 | Submitting author ORCID | MET | `0000-0002-4651-9890` for C. Haden Smith |
| 28 | Author affiliations | MET | 3 affiliations listed |
| 29 | All references have DOIs where available | MET | 21 DOIs + 7 URLs for web-only resources |
| 30 | Code example in paper | MET | Lines 85-113 (bootstrap uncertainty analysis) |

### Repository Metadata

| # | Requirement | Status | Evidence |
|---|-------------|--------|----------|
| 31 | CITATION.cff | MET | `CITATION.cff` (CFF v1.2.0, BSD-3-Clause) |
| 32 | codemeta.json | MET | `codemeta.json` created |
| 33 | Version tag matching submission | REMAINING | Need to create `v2.0.0` tag after merge to `main` |
| 34 | Archived release with DOI | REMAINING | Need Zenodo archive |

---

## Section B: BES Paper Review Comments (Addressed)

| # | Comment | Status | Action |
|---|---------|--------|--------|
| 1 | Co-author ORCIDs (Fields, Gonzalez, Niblett, Beam) | PENDING | Follow up with co-authors separately |
| 2 | "Differential Evolution MCMC" → "Adaptive Differential Evolution MCMC" (line 53) | DONE | Clarifies DE-MCz vs DE-MC |
| 3 | "needed for reliable calibration" → "commonly applied for Bayesian inference" (line 61) | DONE | More accurate characterization |

---

## Section C: Pre-Release Verification

- [ ] Run full test suite across all target frameworks:
  ```
  dotnet test --framework net8.0
  dotnet test --framework net9.0
  dotnet test --framework net10.0
  dotnet test --framework net481
  ```
- [ ] Confirm all 1,006+ tests pass on every target
- [ ] Build the NuGet package locally:
  ```
  dotnet pack Numerics/Numerics.csproj -c Release
  ```
- [ ] Verify the .nupkg contains assemblies for all 4 TFMs

---

## Section D: Release Roadmap (Step-by-Step)

### Step 1: Create Pull Request to `main`
- [ ] Push latest `bugfixes-and-enhancements` to origin
- [ ] Create PR:
  ```
  gh pr create --base main --head bugfixes-and-enhancements \
    --title "v2.0.0: Major update" --body-file RELEASE_NOTES.md
  ```
- [ ] Wait for CI (Integration.yml) to pass
- [ ] Review and merge PR

### Step 2: Tag and GitHub Release
- [ ] After merge:
  ```
  git checkout main && git pull
  git tag v2.0.0
  git push origin v2.0.0
  ```
- [ ] Create GitHub Release:
  ```
  gh release create v2.0.0 --title "v2.0.0 — Major Update" --notes-file RELEASE_NOTES.md
  ```
- [ ] `Release.yml` auto-triggers → pushes to internal USACE Nexus
- [ ] `NuGetPublish.yml` auto-triggers → pushes to nuget.org

### Step 3: One-Time NuGet.org Setup (before first release)
- [ ] Create/verify nuget.org account at https://www.nuget.org/
- [ ] Generate API key at https://www.nuget.org/account/apikeys:
  - Name: `GitHub Actions - Numerics`
  - Expiration: 365 days
  - Glob pattern: `RMC.Numerics`
  - Scopes: "Push new packages and package versions"
- [ ] Add secret to GitHub repo at https://github.com/USACE-RMC/Numerics/settings/secrets/actions:
  - Name: `NUGET_ORG_API_KEY`
  - Value: the API key from above

### Step 4: Verify NuGet Package
- [ ] Check https://www.nuget.org/packages/RMC.Numerics/ (may take 10-15 min to index)
- [ ] Test installation:
  ```
  dotnet new console -o TestInstall && cd TestInstall
  dotnet add package RMC.Numerics --version 2.0.0
  dotnet build
  ```

### Step 5: Zenodo Archival
- [ ] Go to https://zenodo.org and log in with GitHub
- [ ] Enable the `USACE-RMC/Numerics` repository in Zenodo's GitHub integration
- [ ] Zenodo will automatically archive the GitHub Release and mint a DOI
- [ ] Copy the Zenodo DOI badge and add it to `README.md`
- [ ] Update `CITATION.cff` with the Zenodo DOI if desired

### Step 6: Submit to JOSS
- [ ] Go to https://joss.theoj.org/papers/new
- [ ] Enter the repository URL: `https://github.com/USACE-RMC/Numerics`
- [ ] Enter the Zenodo archive DOI
- [ ] Confirm software version matches the tagged release
- [ ] Submit the paper

### Step 7: Post-Submission
- [ ] Add the JOSS status badge to `README.md` once the review issue is created:
  ```markdown
  [![status](https://joss.theoj.org/papers/<DOI>/status.svg)](https://joss.theoj.org/papers/<DOI>)
  ```
- [ ] Monitor the review issue for editor/reviewer feedback
- [ ] Address any review comments promptly (JOSS reviews typically take 4-8 weeks)

---

## Section E: v2.0.0 Release Notes

### v2.0.0 — Major Update

This is a major update to Numerics with 274 files changed, 24,476 insertions, and 4,400 deletions since v1.0.0. Highlights include new distributions, improved MCMC inference, enhanced numerical methods, and comprehensive documentation.

#### New Distributions
- Dirichlet distribution (multivariate)
- Multinomial distribution (multivariate)
- Multivariate Student-t distribution
- Student-t copula

#### Bayesian Inference & MCMC
- Improved Gelman-Rubin convergence diagnostics
- Refactored Noncentral-T to use Brent.Solve
- Enhanced MCMC sampler reliability and convergence

#### Numerical Methods
- Linear algebra enhancements
- Root-finding improvements
- ODE solver improvements
- Improved adaptive integration

#### Data & Statistics
- Time series download improvements
- Hypothesis test enhancements
- Enhanced parameter estimation methods
- Autocorrelation and convergence diagnostics improvements

#### Optimization
- Comprehensive correctness improvements across all optimizers

#### Machine Learning
- Documentation consolidation
- Code quality improvements

#### Infrastructure & Documentation
- Added .NET 10.0 target framework (now targets net8.0, net9.0, net10.0, net481)
- Zero runtime dependencies maintained
- 1,006+ unit tests validated against published references
- JOSS paper and metadata (CITATION.cff, codemeta.json)
- CONTRIBUTING.md and CODE_OF_CONDUCT.md
- 25+ documentation files covering all library capabilities
- NuGet publishing workflow for nuget.org

---

## Notes

- **Word limit**: JOSS papers should be 750-1,750 words (https://joss.readthedocs.io/en/latest/paper.html). Current paper is ~1,215 words.
- **"Jery R. Stedinger"**: This spelling in `paper.bib` is correct (confirmed on USGS publications). Be prepared to explain if a reviewer questions it.
- **England 2019 vs 2018**: Bulletin 17C was originally published March 2018; the `2019` date refers to the v1.1 revision. Both are acceptable in the literature.
- **ter Braak citation**: The 2008 paper describes DE-MCzs (with snooker updater). The library has both `DEMCz` and `DEMCzs` classes. The citation is appropriate since the 2008 paper supersedes the 2006 original.
- **NuGet API key expiration**: The nuget.org API key expires after 365 days max. Set a calendar reminder to regenerate it.
- **Versioning**: Version is derived from git tags (e.g., `v2.1.0` → package version `2.1.0`). The `.csproj` version is a fallback for local builds only.
