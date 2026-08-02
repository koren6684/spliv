# Sensitivity-path optimization audit

Audit date: 2026-08-02. This baseline was recorded before implementation
changes on branch `optimize-sensitivity-path` at commit
`c03d73739eaa814e423bb7384b9befb6f619f202`. The working tree was clean and
`DESCRIPTION` reported version 0.2.1.

## Baseline tests

`devtools::test()` completed with 227 passes, zero failures, zero warnings,
and zero skips in 2.8 seconds.

## Current call chain

1. `spliv_sensitivity_path()` validates the outer delta grid and `...`, then
   calls `.fit_for_sensitivity_delta()` once at delta zero and once for every
   nonzero outer-grid value. A zero already present in the path reuses the
   separately fitted baseline.
2. `.fit_for_sensitivity_delta()` calls the exported `spliv()` from scratch.
3. `spliv()` resolves the method and delegates to `.spliv_impl()`.
4. `.spliv_impl()` calls `.iv_parse()`. That helper parses the two-part IV
   formula, checks referenced variables, constructs the complete-case mask,
   subsets the data, and builds outcome, treatment/control, and instrument
   matrices with `model.matrix()`.
5. `.spliv_impl()` derives the common exogenous-control matrix. With fixed
   effects it removes constants, constructs the fixed-effect frame, and calls
   `.demean_fixest()` or `.demean_lfe()` on the outcome and all model matrices.
6. For clustered inference, `.get_cluster_id()` evaluates the cluster formula
   on the already filtered data (or validates a pre-aligned vector).
7. `.instrument_residual_sds()` partials each excluded instrument against the
   residualized common controls and obtains the scale used by
   `scale_instrument = "residual_sd"`.
8. For patterned sensitivity, `.prepare_violation_pattern()` evaluates and
   normalizes the pattern on the estimation sample, identifies the single
   treatment/instrument, applies instrument scaling, and constructs the
   invariant direct-effect regressor.
9. UCI constructs the inner gamma/theta grid and calls `.sp_uci_mats()`. For
   each inner point, that helper forms the adjusted outcome and calls
   `.iv_2sls_mats()`.
10. LTZ constructs the delta-specific prior and calls `.sp_ltz_mats()`, which
    calls `.iv_2sls_mats()` for the unchanged baseline IV fit and then combines
    it with the direct-effect prior.
11. `.iv_2sls_mats()` recomputes `Z'Z`, its inverse, `X'Z`, `X'PzX`, the
    bread, residuals, and the IID/HC1/cluster covariance ingredients.
12. The path extracts identically structured rows, computes zero containment,
    and derives grid-based tipping points after all outer-grid fits finish.

## Repeated invariant work identified

For a path with `D` nonzero values, the baseline-plus-loop structure repeats
formula parsing, variable checking, complete-case filtering, model-matrix
construction, fixed-effect-frame construction, fixed-effect residualization,
cluster alignment, residual instrument scaling, and pattern evaluation
`D + 1` times. If the supplied grid omits zero, the extra baseline fit is still
performed for path metadata and tipping-point comparisons.

The matrix layer also repeats invariant cross-products and matrix inversions.
UCI repeats them for every point of every inner gamma/theta grid. LTZ repeats
the identical baseline IV fit, covariance, and direct-effect loading matrix for
every outer delta; only its prior covariance varies in the standard path.

The adjusted outcome itself is already formed after fixed-effect
residualization. The optimization must retain the package's current
direct-effect construction and use the linear identity
`M_A(Y - theta G) = M_A Y - theta M_A G` to keep the invariant residualized
components outside the outer loop without changing the estimand.

## Planned behavior-preserving boundary

The optimization will introduce only private prepared objects. The public
signature, returned columns/classes/attributes, delta and inner-grid meanings,
pattern normalization, instrument scaling, covariance definitions, warnings,
and errors will remain unchanged. Confirmatory and exploratory BPE are outside
the path and will not be altered.

## Implemented strategy

The final implementation adds a private `spliv_prepared_design` that stores the
single parsed/filtered estimation sample, residualized outcome and model
matrices, aligned cluster id, residual instrument scale, evaluated pattern,
and warnings. A private `spliv_prepared_iv_design` stores the invariant 2SLS
cross-products, inverse, bread, and cluster grouping.

`spliv_sensitivity_path()` prepares these objects once and passes them to each
delta-specific internal fit. UCI reuses the same 2SLS design across outer and
inner grids while continuing to recompute the delta/theta-specific adjusted
outcome, residuals, and covariance meat. LTZ additionally reuses the unchanged
baseline IV fit and direct-effect loading matrix; only its delta-specific prior
combination is recalculated. The exported function signature and returned path
structure are unchanged.

Test-only instrumentation records actual preparation events. On a clustered,
two-way-FE patterned path with nine deltas, formula/model preparation,
complete-case filtering, fixed-effect residualization, cluster alignment,
instrument scaling, pattern evaluation, and invariant IV-design preparation
each occurred exactly once.

## Final numerical evidence

A transparent reference helper retains the former baseline-plus-repeated-
`spliv()` behavior without passing a prepared design. Comparisons cover no,
one-way, and two-way fixed effects; IID, HC1, and clustered covariance; uniform
and patterned UCI/LTZ; missing outcome/treatment/instrument values; missing
pattern error behavior; cluster alignment; formula/function/vector/column
patterns; raw and residual-SD scaling; instruments multiplied by 100 and 0.01;
short/long/zero-containing grids; classes, columns, metadata, zero containment,
and tipping points.

The maximum finite numerical difference over estimates, standard errors,
interval endpoints, baselines, and tipping-point inputs was exactly zero. The
final suite contains 1,756 passing expectations: 227 pre-existing and 1,529
new, with zero failures, warnings, or skips.

## Final validation

- `roxygen2::roxygenise()`: PASS; no generated public API change.
- `devtools::check(document = FALSE)`: PASS, 0 errors, 0 warnings, 0 notes.
- `R CMD build .`: PASS for `spliv_0.2.1.tar.gz`.
- Exact-tarball `R CMD check --as-cran`: PASS under the installed
  `en_US.UTF-8` locale with 0 errors, 0 warnings, and 1 incoming-feasibility
  NOTE (`Days since last update: 4`). Tests, examples, vignette rebuild, PDF
  manual, and HTML manual passed.
- The first exact check under the desktop process's unsupported inherited
  `C.UTF-8` locale ended with 1 error, 1 warning, and 1 note because R emitted
  locale startup warnings during DESCRIPTION checking. Its check directory was
  preserved under `/tmp`; the identical tarball then passed under the valid
  installed locale.

`DESCRIPTION` remains version 0.2.1. `NAMESPACE` and all exported signatures
are unchanged. No commit, tag, push, release, CRAN submission, benchmark,
simulation, replication, article build, or change outside the package
repository was performed.

## Remaining performance work

The outer-loop parsing, filtering, demeaning, scaling, and design setup are no
longer bottlenecks. UCI must still form each adjusted outcome and compute its
coefficient, residuals, and IID/HC1/cluster covariance contribution at every
unique inner-grid point. It also still performs a small linear solve for every
point even though the bread is cached, and the path still constructs a full
`spliv_fit` list before extracting each row. Those are candidates for the
separate benchmark-driven pass; none was changed speculatively here.

The package is ready for a benchmark-only rerun.
