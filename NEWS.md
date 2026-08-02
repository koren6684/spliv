# spliv 0.2.1

- `spliv_sensitivity_path()` now reuses invariant formula parsing,
  complete-case filtering, model matrices, fixed-effect residualization,
  cluster alignment, pattern evaluation, instrument scaling, and 2SLS design
  calculations across its delta grid. This is an internal, numerically
  equivalent optimization: the public API and statistical results are not
  intended to change.
- `spliv()` is the sole recommended primary estimator.
- Fit objects now use class `spliv_fit`.
- Accidental legacy aliases present in the 0.1.x interface were removed from
  the public API, including the old `plausexog()`/`conley_*()` names and
  implementation helper exports.
- `spliv()` now defaults to conventional UCI at `delta = 0` when no sensitivity
  method or bound is supplied, and obsolete BPE arguments were removed.
- Confirmatory BPE equivalence intervals are now compared on the documented,
  unit-invariant scale: residual treatment standard deviations per one-
  residual-SD instrument shift. Raw and standardized diagnostics are retained
  separately. This corrects 0.1.x eligibility calculations that compared a raw
  coefficient interval with a scale-aware margin.
- Confirmatory BPE now requires a non-empty transportability rationale and uses
  the aligned `sampling` or `conservative` covariance modes in standalone
  validation and final estimation.
- Exploratory subset diagnostics retain first-stage F-statistics as descriptive
  information but no longer expose an F-threshold selection control.
- Numerical-equivalence, BPE invariance, fixed-effect, clustered-inference,
  edge-case, and plotting tests were expanded. No treatment-effect formula or
  sensitivity estimand was intentionally changed apart from the BPE
  equivalence-scale correction above.
- A finalized 72-case equivalent-computation benchmark was added under
  `dev/`. Corresponding sensitivity-interval endpoints agreed to numerical
  precision. The optimized sensitivity-path implementation was faster than
  repeated adjusted-outcome refitting in every benchmarked design. Allocated
  memory was substantially lower in many cases but was not uniformly lower,
  including one large design with approximately 25% greater allocation.

# spliv 0.2.0

- GitHub-only public API cleanup release. It was superseded by version 0.2.1
  before the corresponding CRAN update.

# spliv 0.1.1

- Include generated vignette metadata in source-package builds.
- Complete initial CRAN-readiness and portability checks.

# spliv 0.1.0

Initial public staging release. This version provides uniform and patterned
UCI/LTZ sensitivity analysis, sensitivity paths and tipping points, fixed-effect
residualization, confirmatory BPE design/validation, and compatibility wrappers
for the earlier plausibly-exogenous API.
