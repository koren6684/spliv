# Patterned Sensitivity Analysis for Plausibly Exogenous IV

Main estimator for patterned sensitivity analysis of exclusion
violations in fixed-effect or residualized IV designs.

## Usage

``` r
plausexog(
  formula,
  data,
  fe = NULL,
  fe_engine = c("fixest", "lfe"),
  vcov = c("iid", "hc1", "cluster"),
  cluster = NULL,
  method = c("ltz", "uci", "bpe"),
  prior = NULL,
  delta = NULL,
  violation_pattern = NULL,
  bpe = FALSE,
  bpe_design = NULL,
  bpe_spec = list(design = NULL, subset = NULL, z_names = NULL),
  bpe_kappa = 1,
  bpe_omega = NULL,
  bpe_min_n_S = 2000,
  bpe_min_clusters_S = 30,
  bpe_max_F_S = NULL,
  bpe_min_varZ_S = 1e-06,
  bpe_equiv_margin = NULL,
  bpe_equiv_level = 0.95,
  bpe_transport = c("sampling", "none", "conservative"),
  bpe_transport_kappa = 0,
  bpe_not_applicable = c("na", "error"),
  scale_instrument = c("residual_sd", "none"),
  grid = list(),
  ...
)

plausexog_iv(...)

spliv(
  formula,
  data,
  fe = NULL,
  fe_engine = c("fixest", "lfe"),
  vcov = c("iid", "hc1", "cluster"),
  cluster = NULL,
  method = c("ltz", "uci", "bpe"),
  prior = NULL,
  delta = NULL,
  violation_pattern = NULL,
  bpe = FALSE,
  bpe_design = NULL,
  bpe_spec = list(design = NULL, subset = NULL, z_names = NULL),
  bpe_kappa = 1,
  bpe_omega = NULL,
  bpe_min_n_S = 2000,
  bpe_min_clusters_S = 30,
  bpe_max_F_S = NULL,
  bpe_min_varZ_S = 1e-06,
  bpe_equiv_margin = NULL,
  bpe_equiv_level = 0.95,
  bpe_transport = c("sampling", "none", "conservative"),
  bpe_transport_kappa = 0,
  bpe_not_applicable = c("na", "error"),
  scale_instrument = c("residual_sd", "none"),
  grid = list(),
  ...
)
```

## Arguments

- formula:

  IV formula `y ~ X | Z`.

- data:

  Data frame.

- fe:

  Optional one-sided formula of fixed effects.

- fe_engine:

  FE demeaning engine, one of `"fixest"` or `"lfe"`.

- vcov:

  One of `"iid"`, `"hc1"`, or `"cluster"`.

- cluster:

  Cluster ids or one-sided formula, required when `vcov = "cluster"`.

- method:

  One of `"ltz"`, `"uci"`, or `"bpe"`.

- prior:

  Optional prior list with `mu` and `Omega` (or `omega`) for LTZ. When
  `violation_pattern` is supplied, patterned LTZ currently requires a
  scalar prior over the pattern coefficient.

- delta:

  Optional non-negative scalar sensitivity magnitude. With
  `method = "uci"`, `delta` implies theta bounds `[-delta, +delta]`
  unless explicit bounds are supplied in `grid`. With `method = "ltz"`
  and no explicit `prior`, `delta` induces a zero-mean normal LTZ prior.

- violation_pattern:

  Optional
  [`spliv_pattern()`](https://koren6684.github.io/spliv/reference/spliv_pattern.md)
  object describing how the direct effect of the instrument may vary
  across observations. If omitted, LTZ/UCI retain the package's
  backward-compatible uniform direct-effect behavior. This argument is
  currently supported for LTZ and UCI, but not for confirmatory BPE.

- bpe:

  Logical; when `TRUE`, learn prior moments from a confirmatory
  [`bpe_design()`](https://koren6684.github.io/spliv/reference/bpe_design.md)
  and run LTZ.

- bpe_design:

  Optional
  [`bpe_design()`](https://koren6684.github.io/spliv/reference/bpe_design.md)
  object for confirmatory BPE.

- bpe_spec:

  Optional list. Prefer `bpe_spec = list(design = my_design)`. For
  backward compatibility, `bpe_spec$subset` may also be supplied, but it
  must represent an explicit researcher-supplied subset and should be
  paired with a non-empty `rationale`, explicit `pre_specified = TRUE`,
  and any relevant metadata such as `variables_used` or `subset_type`.
  Supply exactly one confirmatory subset source: `bpe_design`,
  `bpe_spec$design`, or `bpe_spec$subset`. Exploratory `subset_rule`
  inputs are not accepted for confirmatory estimation.

- bpe_kappa:

  Positive scalar multiplier applied to the confirmatory BPE covariance
  after transport adjustment.

- bpe_omega:

  Deprecated. Confirmatory BPE now uses the full reduced-form covariance
  from the subset together with `bpe_transport`.

- bpe_min_n_S:

  Minimum subset size required for BPE eligibility. Default `2000`.

- bpe_min_clusters_S:

  Minimum number of clusters required in subset `S` when
  `vcov = "cluster"`. Default `30`.

- bpe_max_F_S:

  Deprecated. The first-stage F-statistic is reported for diagnostics
  only and no longer determines confirmatory BPE eligibility.

- bpe_min_varZ_S:

  Minimum residualized instrument variance required in subset `S`.
  Default `1e-6`.

- bpe_equiv_margin:

  Researcher-specified first-stage equivalence margin. Confirmatory BPE
  currently supports exactly one instrument.

- bpe_equiv_level:

  Confidence level for the first-stage equivalence check.

- bpe_transport:

  One of `"none"`, `"sampling"`, or `"conservative"`.

- bpe_transport_kappa:

  Non-negative scalar controlling the conservative transport covariance
  inflation.

- bpe_not_applicable:

  Behavior when subset diagnostics fail. One of `"na"` (default) to
  return NA estimates, or `"error"` to stop.

- scale_instrument:

  One of `"residual_sd"` (default) or `"none"`.

- grid:

  List controlling UCI bounds or other tuning parameters. For
  backward-compatible scalar UCI, if `grid$delta` is supplied and
  `grid$gmin`/`grid$gmax` are omitted, the package interprets `delta` as
  a direct-effect bound of `[-delta, +delta]` under the chosen
  `scale_instrument`. When `violation_pattern` is supplied, `grid$delta`
  instead refers to theta bounds over the pattern-scaled direct effect.

- ...:

  Reserved.

## Value

Object of class `plausexog_fit`.

## Details

`spliv` implements patterned sensitivity analysis for exclusion
violations in IV designs with fixed effects or other residualization
steps. Researchers can supply a theoretically motivated
[`spliv_pattern()`](https://koren6684.github.io/spliv/reference/spliv_pattern.md)
object to specify where direct effects of an instrument are expected to
be larger or smaller, and the package then scales LTZ/UCI sensitivity
along that pattern.

The package does not estimate an unrestricted direct-effect field.
Instead, researchers supply a structured pattern and ask whether
conclusions survive direct effects scaled along that pattern.

In applied work, users should usually vary `delta` over a range with
[`spliv_sensitivity_path()`](https://koren6684.github.io/spliv/reference/spliv_sensitivity_path.md)
rather than report one arbitrary sensitivity value.

Patterned sensitivity currently supports one endogenous treatment, one
excluded instrument, and one researcher-specified pattern at a time.

Confirmatory BPE is not a subgroup-search procedure. The researcher must
supply a pre-specified
[`bpe_design()`](https://koren6684.github.io/spliv/reference/bpe_design.md)
object, the package validates that subset, and BPE proceeds only if the
confirmatory eligibility checks pass.

The first-stage F-statistic is still reported for diagnostics, but
confirmatory BPE eligibility is determined by the pre-specification
checks, subset size, cluster count, residualized instrument variation,
and a first-stage equivalence interval.

BPE reduced-form covariance is propagated as a full covariance matrix
and can optionally be inflated via `bpe_transport`.

## Deprecated Wrappers

`plausexog()` and `plausexog_iv()` are deprecated compatibility aliases
for `spliv()`.

## Examples

``` r
set.seed(1)
d <- data.frame(y = rnorm(60), x = rnorm(60), z = rnorm(60), w = rnorm(60))
fit <- spliv(y ~ x + w | z + w, d, method = "uci", delta = 0.1)
fit$estimates
#>          term    conf.low conf.high
#> 1 (Intercept)  -0.8544988  1.363666
#> 2           x -10.1507140  7.661939
#> 3           w  -0.5514807  0.410585
```
