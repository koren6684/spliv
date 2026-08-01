# Final SPLIV 0.2.0 R Journal package audit

Audit date: 2026-07-31. Status labels refer to the current uncommitted staging
tree. `PASS` means the requested local evidence is present; `BLOCKED` means a
maintainer or long-running production action is still required; `FAIL` means a
completed validation did not meet its criterion.

| Criterion | Status | Evidence |
|---|---|---|
| CRAN availability of stable package | PASS | The official CRAN page, <https://cran.r-project.org/package=spliv>, reported stable version 0.1.1 when checked on 2026-07-31. The staged 0.2.0 has not been submitted. |
| Intended 0.2.0 API | PASS | `NAMESPACE` and a namespace test contain exactly nine canonical and seven advanced exports. No removed alias is exported or documented. |
| Package substance | PASS | The package contains separate parsing, matrix, FE, BPE, estimator, and plotting modules; 16 documented exports; a workflow vignette; and a substantive test suite. |
| Comparison with related packages | BLOCKED | Numerical comparisons with `fixest` establish baseline-IV equivalence and a fair adjusted-outcome benchmark. A manuscript-level discussion comparing SPLIV's scope with related sensitivity-analysis packages has not yet been written or peer reviewed. |
| Unit tests | PASS | Final source tests: 227 passed, zero failures, warnings, or skips. |
| Test count and coverage areas | PASS | Tests cover baseline IV, UCI, LTZ, patterns, BPE, FE backends, covariance choices, scale invariance, transport, plotting, public API, and errors/edge cases. A formal line-coverage percentage was not required or computed. |
| Documentation for every export | PASS | Every approved export has a title, description, usage, arguments, value, scope details where needed, and generated Rd help. |
| Examples for every export | PASS | All 16 exports have small offline examples; exact-tarball `R CMD check` ran all examples successfully. |
| Full workflow vignette | PASS | `vignettes/getting-started.Rmd` uses a coherent patterned synthetic DGP and demonstrates baseline IV, UCI, patterned UCI/LTZ, paths, tipping points, and confirmatory BPE. It rebuilt during check. |
| Version control | PASS | The package is in a public Git repository on branch `main`; the baseline commit is `41993897f83f9388abd2cd13368c2e36ae91aefa`. This audit intentionally leaves a dirty, uncommitted release-candidate tree. |
| Issue tracker | PASS | <https://github.com/koren6684/spliv/issues> returned HTTP 200 on 2026-07-31. |
| Dependency discipline | PASS | Imports remain `stats`, `graphics`, and `fixest`; optional FE, plotting, testing, and vignette packages are in Suggests. Dependency checks passed. |
| Coherent naming | PASS | The recommended interface is consistently `spliv*`/`bpe_*`; genuinely distinct low-level functions use `sp_*`. Historical aliases are absent from the current interface. |
| S3 methods | PASS | New fits have class exactly `spliv_fit`; print/plot methods and the sensitivity-path plot method are registered. Unexported old-class methods remain only for serialized 0.1.x objects. |
| Error and warning behavior | PASS | Unsupported designs, singular matrices, missing variables, invalid patterns, absent plot terms, and ineligible BPE designs produce tested informative outcomes. Exploratory BPE warns that it is non-confirmatory. |
| No unnecessary file I/O | PASS | Estimators do not write files or access the network. Examples and vignette are offline; plotting writes only when the caller opens a device. |
| pkgdown website | PASS | A clean local site build completed, with approved reference sections and no obsolete reference pages. The public URL returned HTTP 200; deployment of the uncommitted changes was not attempted. |
| Numerical equivalence | PASS | SPLIV coefficients and IID/HC1/cluster SEs match equivalent `fixest::feols` IV fits with no FE and two-way FE to `1e-9`; manual adjusted-outcome UCI endpoints match to `1e-9`. |
| BPE scale invariance | PASS | Tests replacing `z` by `100*z` and `0.01*z` preserve standardized first-stage effects/CIs, equivalence decisions, eligibility, and final BPE intervals. |
| BPE transportability logic | PASS | A nonempty transportability rationale is required. `sampling` transports the estimated sampling covariance; `conservative` multiplies it by `1 + kappa`. Standalone validation and final fitting use identical priors. |
| Fixed-effect engine consistency | PASS | Guarded tests compare `fixest` and `lfe` two-way-FE results at `1e-10`; the installed local environment ran the test without a skip. |
| Clustered inference | PASS | Baseline clustered IV is compared directly with `fixest`; character and formula cluster inputs agree, missing cluster rows are filtered, and insufficient BPE clusters are tested. |
| Edge-case handling | PASS | Tests cover missingness, weak instruments, collinearity, zero instrument variation, singular models, invalid/all-zero/NA patterns, multiple treatments/patterns, small subsets, and plot-term errors. |
| Performance and scaling | PASS | The 72-case full benchmark completed with zero errors and endpoint differences at most `5.24e-14`. SPLIV is faster in no-FE and smaller FE cases but slower/more allocation-intensive for every `n=50,000` two-way-FE case; no blanket speed claim remains. |
| Source tarball cleanliness | PASS | The exact 0.2.0 tarball excludes `dev/`, `docs/`, Git metadata, local/restricted data, logs, checks, benchmark output, and release PDFs. Standard built vignette files are present. |
| Reference manual | PASS | PDF and HTML manuals built. All 20 PDF pages were rendered and visually inspected; signatures, arguments, examples, and index contain the approved API without stale aliases. |
| Reproducibility repository | BLOCKED | Clean-room installed-package environment, smoke, and all pilot families passed with zero failed tasks. The full 1,000-replicate simulation references have not been regenerated after the BPE scale correction, and `renv.lock` points to an older tagged commit. |
| Koren article use case | PASS | The authorized external data passed structural/checksum checks. The full public runner completed all 21 delta values for both crops and patterns, reproduced the expected UCI tipping points (`0.01` uniform and `0.03` sparse/bare), and produced the corrected confirmatory BPE result: maize eligible, wheat ineligible. The compact R Journal table/figure builder then completed and its six CSV and two figure references match byte-for-byte. No empirical data were copied into either repository. |
| Public data availability | PASS | Koren's permanent Harvard Dataverse archive is DOI `10.7910/DVN/Q3UISS`; released metadata identifies CC0 1.0 and an unrestricted archive. Data remain external and checksum-gated. Lelkes is optional and its redistribution terms remain unconfirmed. |
| Article runtime feasibility | PASS | The full Koren run completed locally in approximately 57 minutes 38 seconds (11:08:23--12:06:01), including both crops, two patterns, and 21-value UCI/LTZ grids. The subsequent compact article-output build completed successfully. This establishes practical local runtime for the principal use case; an end-to-end R Journal manuscript build remains a submission blocker. |
| Remaining work before CRAN 0.2.0 | BLOCKED | Resolve the already-existing `v0.2.0` Git tag, commit the audited tree, run hosted multi-platform CI, review the expected recent-update NOTE, and rerun the exact final tarball check after any release metadata change. |
| Remaining work before R Journal submission | BLOCKED | Complete/freeze the corrected full simulation results; regenerate `renv.lock` to the actual immutable release; finish the related-software comparison and manuscript; run its end-to-end build; and archive the exact article artifacts. |

## Package validation details

- `roxygen2::roxygenise()`: PASS.
- `devtools::test()`: 227 passed, 0 failed, 0 warnings, 0 skips.
- `pkgdown::build_site()`: PASS; Pandoc reported only its upstream
  `--highlight-style` deprecation message.
- `R CMD build .`: PASS.
- Exact-tarball `R CMD check --as-cran`: 0 errors, 0 warnings, 1 NOTE
  (maintainer/recent update; two days since last update). Examples, tests,
  vignette rebuild, PDF manual, and HTML manual passed.
- An initial check under this desktop session's invalid `C.UTF-8` environment
  ended with 1 error, 1 warning, and 1 note because R emitted locale warnings
  during metadata checking. The same tarball passed under the installed
  `en_US.UTF-8` locale; the failed attempt is not concealed as a package issue.

## Scientific change record

The only intended numerical correction is the confirmatory BPE equivalence
scale. With residual-SD scaling, eligibility now compares the first-stage CI
after multiplying raw coefficients by `sd(z_residual)/sd(x_residual)` with the
dimensionless equivalence margin. UCI, LTZ, and treatment-effect formulas were
not changed. The deterministic pilot's BPE eligibility/fit-success rates
changed from `0.1/0.1/0.0` to `0.0/0.1/0.0`; patterned and subgroup-search
pilot summaries were unchanged. Koren and Lelkes exact changes are recorded in
the release audit.

## Release-reference blocker

The remote lightweight `v0.2.0` tag already resolves to baseline commit
`41993897f83f9388abd2cd13368c2e36ae91aefa`, not this final uncommitted tree.
The checked-in `renv.lock` resolves that same commit. No tag or lockfile was
changed in this pass. The maintainer must choose an immutable release/version
strategy and regenerate the lockfile; silently treating the current tag as the
audited release would be scientifically incorrect.
