# spliv 0.2.0

- `spliv()` is the sole recommended primary estimator.
- Fit objects now use class `spliv_fit`.
- Accidental legacy aliases introduced in 0.1.1 were removed from the public
  API, including the old `plausexog()`/`conley_*()` names and implementation
  helper exports.
- No statistical estimators or numerical results were intentionally changed.

# spliv 0.1.1

- Include generated vignette metadata in source-package builds.
- Complete initial CRAN-readiness and portability checks.

# spliv 0.1.0

Initial public staging release. This version provides uniform and patterned
UCI/LTZ sensitivity analysis, sensitivity paths and tipping points, fixed-effect
residualization, confirmatory BPE design/validation, and compatibility wrappers
for the earlier plausibly-exogenous API.
