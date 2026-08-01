# Patterned Sensitivity Analysis for Plausibly Exogenous IV

Main estimator for patterned sensitivity analysis of exclusion
violations in fixed-effect or residualized IV designs.

## Usage

``` r
spliv(
  formula,
  data,
  fe = NULL,
  fe_engine = c("fixest", "lfe"),
  vcov = c("iid", "hc1", "cluster"),
  cluster = NULL,
  method = c("uci", "ltz", "bpe"),
  prior = NULL,
  delta = NULL,
  violation_pattern = NULL,
  bpe_design = NULL,
  bpe_kappa = 1,
  bpe_min_n_S = 2000,
  bpe_min_clusters_S = 30,
  bpe_min_varZ_S = 1e-06,
  bpe_equiv_margin = NULL,
  bpe_equiv_level = 0.95,
  bpe_transport = c("sampling", "conservative"),
  bpe_transport_kappa = 0,
  bpe_not_applicable = c("na", "error"),
  scale_instrument = c("residual_sd", "none"),
  grid = list()
)
```

## Arguments

- formula:

  IV formula `y ~ X | Z`. Ordinary exogenous controls must appear on
  both sides, for example `y ~ x + w | z + w`.

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

  One of `"uci"`, `"ltz"`, or `"bpe"`. The default is the conventional
  UCI analysis.

- prior:

  Optional prior list with `mu` and `Omega` (or `omega`) for LTZ. When
  `violation_pattern` is supplied, patterned LTZ currently requires a
  scalar prior over the pattern coefficient.

- delta:

  Optional non-negative scalar sensitivity magnitude. With
  `scale_instrument = "residual_sd"`, delta is the direct outcome effect
  of a one-residual-SD instrument shift. With
  `scale_instrument = "none"`, it is the raw direct-effect coefficient.
  For UCI, delta supplies symmetric bounds `[-delta, +delta]`; for LTZ
  without an explicit prior, it is the standard deviation of a zero-mean
  normal direct-effect prior.

- violation_pattern:

  Optional
  [`spliv_pattern()`](https://koren6684.github.io/spliv/reference/spliv_pattern.md)
  object describing how the direct effect of the instrument may vary
  across observations. If omitted, LTZ/UCI use the uniform direct-effect
  behavior. This argument is currently supported for LTZ and UCI, but
  not for confirmatory BPE.

- bpe_design:

  A pre-specified
  [`bpe_design()`](https://koren6684.github.io/spliv/reference/bpe_design.md)
  object for confirmatory BPE. A non-empty subset rationale and
  transportability rationale are required when `method = "bpe"`.

- bpe_kappa:

  Positive scalar multiplier applied to the confirmatory BPE covariance
  after transport adjustment.

- bpe_min_n_S:

  Minimum subset size required for BPE eligibility. Default `2000`.

- bpe_min_clusters_S:

  Minimum number of clusters required in subset `S` when
  `vcov = "cluster"`. Default `30`.

- bpe_min_varZ_S:

  Minimum residualized instrument variance required in subset `S`.
  Default `1e-6`.

- bpe_equiv_margin:

  Researcher-specified first-stage equivalence margin. With
  `scale_instrument = "residual_sd"`, the margin is measured in residual
  treatment standard deviations per one-residual-SD instrument shift.
  With `scale_instrument = "none"`, it is on the raw first-stage
  coefficient scale. Confirmatory BPE currently supports one treatment
  and one instrument.

- bpe_equiv_level:

  Confidence level for the first-stage equivalence check.

- bpe_transport:

  One of `"sampling"` or `"conservative"`. Sampling uses the estimated
  reduced-form sampling covariance; conservative adds the inflation
  controlled by `bpe_transport_kappa`.

- bpe_transport_kappa:

  Non-negative scalar controlling the conservative transport covariance
  inflation.

- bpe_not_applicable:

  Behavior when subset diagnostics fail. One of `"na"` (default) to
  return NA estimates, or `"error"` to stop.

- scale_instrument:

  One of `"residual_sd"` (default) or `"none"`.

- grid:

  List controlling UCI bounds or other tuning parameters. For scalar
  UCI, if `grid$delta` is supplied and `grid$gmin`/`grid$gmax` are
  omitted, the package interprets `delta` as a direct-effect bound of
  `[-delta, +delta]` under the chosen `scale_instrument`. When
  `violation_pattern` is supplied, `grid$delta` instead refers to theta
  bounds over the pattern-scaled direct effect.

## Value

Object of class `spliv_fit`.

## Details

UCI is the union of conventional IV confidence intervals over the
specified direct-effect bounds. LTZ propagates a local-to-zero prior
mean and covariance for the direct effect into the IV estimate and
uncertainty.

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
object, the package validates that subset, documents why its direct
effect is transportable to the target sample, and BPE proceeds only if
the confirmatory eligibility checks pass.

The first-stage F-statistic is still reported for diagnostics, but
confirmatory BPE eligibility is determined by the pre-specification
checks, subset size, cluster count, residualized instrument variation,
and a first-stage equivalence interval.

BPE reduced-form covariance is propagated as a full covariance matrix
and can optionally be inflated via `bpe_transport`.

## Examples

``` r
set.seed(1)
d <- data.frame(y = rnorm(60), x = rnorm(60), z = rnorm(60), w = rnorm(60))
fit <- spliv(y ~ x + w | z + w, d)
fit$estimates
#>          term   conf.low conf.high
#> 1 (Intercept) -0.3439308 0.6159157
#> 2           x -4.0703086 3.4365025
#> 3           w -0.3025500 0.1376759
```
