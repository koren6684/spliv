# Getting started with spliv

## Baseline IV and uniform sensitivity

The [`spliv()`](https://koren6684.github.io/spliv/reference/spliv.md)
estimator accepts the usual IV formula `y ~ X | Z`. Ordinary exogenous
controls belong on both sides, as in `y ~ x + w | z + w`. With no method
or bound supplied, the baseline is UCI at `delta = 0`; a positive
`delta` allows a bounded uniform direct effect. On the default scale,
`delta` is an outcome-unit direct effect for a one-residual-SD
instrument shift.

``` r

baseline <- spliv(f, d, vcov = "hc1")
baseline$estimates
#>          term   conf.low conf.high
#> 1 (Intercept) -0.1605207 0.1107809
#> 2           x  0.9391974 1.6391314
#> 3           w -0.0851201 0.3415725
uniform <- spliv(f, d, method = "uci", delta = 0.20, vcov = "hc1",
                 grid = list(steps = 9))
uniform$estimates
#>          term   conf.low conf.high
#> 1 (Intercept) -0.2340356 0.1964620
#> 2           x  0.3485695 2.2672576
#> 3           w -0.3722147 0.6281993
```

## A patterned UCI/LTZ analysis

Patterns are explicit objects with a substantive rationale. Here the
direct effect is allowed to be larger where the synthetic exposure is
larger.

``` r

pattern <- spliv_pattern(
  name = "Exposure pattern", pattern = ~ exposure,
  rationale = "The alternative channel is stronger at higher exposure.",
  variables_used = "exposure", pattern_type = "theory_defined",
  normalize = "max_abs"
)
spliv_eval_pattern(pattern, d)[1:5]
#> [1] 0.23328645 0.84214797 0.89724528 0.89549398 0.08382179
patterned_uci <- spliv(f, d, method = "uci", delta = 0.20, vcov = "hc1",
                       violation_pattern = pattern, grid = list(steps = 9))
patterned_ltz <- spliv(f, d, method = "ltz", delta = 0.20, vcov = "hc1",
                       violation_pattern = pattern)
patterned_uci$estimates
#>          term   conf.low conf.high
#> 1 (Intercept) -0.1888420 0.1445298
#> 2           x  0.6828480 1.9134147
#> 3           w -0.2138458 0.4704673
patterned_ltz$estimates
#>          term    estimate  std.error   conf.low conf.high
#> 1 (Intercept) -0.02486991 0.07400209 -0.1699114 0.1201715
#> 2           x  1.28916439 0.30312599  0.6950484 1.8832804
#> 3           w  0.12822619 0.16031139 -0.1859784 0.4424307
```

## Sensitivity paths and tipping points

Paths make the sensitivity grid explicit. A tipping point is reported
only when the supplied interval contains zero at a grid value.

``` r

path <- spliv_sensitivity_path(
  f, d, method = "uci", delta_grid = seq(0, 0.30, by = 0.05),
  vcov = "hc1", violation_pattern = pattern
)
head(path)
#>          term delta method    estimate   conf_low conf_high contains_zero
#> 1 (Intercept)  0.00    uci -0.02486991 -0.1605207 0.1107809          TRUE
#> 2           x  0.00    uci  1.28916439  0.9391974 1.6391314         FALSE
#> 3           w  0.00    uci  0.12822619 -0.0851201 0.3415725          TRUE
#> 4 (Intercept)  0.05    uci -0.02486991 -0.1666758 0.1183381          TRUE
#> 5           x  0.05    uci  1.28916439  0.8790158 1.7040397         FALSE
#> 6           w  0.05    uci  0.12822619 -0.1152324 0.3717291          TRUE
#>       pattern_name   pattern_type violation_pattern_used scale_instrument nobs
#> 1 Exposure pattern theory_defined                   TRUE      residual_sd  240
#> 2 Exposure pattern theory_defined                   TRUE      residual_sd  240
#> 3 Exposure pattern theory_defined                   TRUE      residual_sd  240
#> 4 Exposure pattern theory_defined                   TRUE      residual_sd  240
#> 5 Exposure pattern theory_defined                   TRUE      residual_sd  240
#> 6 Exposure pattern theory_defined                   TRUE      residual_sd  240
#>   se theta_min theta_max baseline_estimate baseline_conf_low baseline_conf_high
#> 1 NA      0.00      0.00       -0.02486991        -0.1605207          0.1107809
#> 2 NA      0.00      0.00        1.28916439         0.9391974          1.6391314
#> 3 NA      0.00      0.00        0.12822619        -0.0851201          0.3415725
#> 4 NA     -0.05      0.05       -0.02486991        -0.1605207          0.1107809
#> 5 NA     -0.05      0.05        1.28916439         0.9391974          1.6391314
#> 6 NA     -0.05      0.05        0.12822619        -0.0851201          0.3415725
#>   crosses_baseline_sign significant_at_level error
#> 1                  TRUE                FALSE  <NA>
#> 2                 FALSE                 TRUE  <NA>
#> 3                  TRUE                FALSE  <NA>
#> 4                  TRUE                FALSE  <NA>
#> 5                 FALSE                 TRUE  <NA>
#> 6                  TRUE                FALSE  <NA>
spliv_tipping_point(path)
#> (Intercept)           x           w 
#>           0          NA           0 
#> attr(,"message")
#>                                             (Intercept) 
#>              "Baseline interval already includes zero." 
#>                                                       x 
#> "No zero crossing occurred on the supplied delta grid." 
#>                                                       w 
#>              "Baseline interval already includes zero."
```

## Confirmatory BPE design and validation

BPE uses an outcome-independent subset specified before examining the
outcome analysis and requires a non-empty transportability rationale.
Validation reports eligibility diagnostics; it does not search across
candidate subgroups. The default `sampling` transport uses the estimated
reduced-form sampling covariance; `conservative` adds a pre-specified
covariance inflation.

``` r

design <- bpe_design(
  name = "Theory-defined inactive subset", subset = ~ inactive,
  rationale = "The treatment channel is absent in the inactive subset.",
  variables_used = "inactive", subset_type = "theory_defined",
  pre_specified = TRUE,
  transportability_rationale = "The subset direct effect is informative for the target sample."
)
# Illustrative synthetic margin: 0.25 residual treatment SD per
# one-residual-SD instrument shift. In substantive work, pre-specify the
# margin rather than tuning it to obtain BPE eligibility.
bpe_margin <- 0.25
validation <- bpe_validate_design(
  f, d, design = design, vcov = "hc1",
  bpe_min_n_S = 40,
  bpe_equiv_margin = bpe_margin
)
validation[c("n_S", "equivalence_passed", "eligibility_passed")]
#> $n_S
#> [1] 120
#> 
#> $equivalence_passed
#> [1] FALSE
#> 
#> $eligibility_passed
#> [1] FALSE

bpe_fit <- spliv(f, d, method = "bpe", bpe_design = design,
                 vcov = "hc1", bpe_min_n_S = 40,
                 bpe_equiv_margin = bpe_margin)
bpe_fit$estimates
#>          term estimate std.error conf.low conf.high
#> 1 (Intercept)       NA        NA       NA        NA
#> 2           x       NA        NA       NA        NA
#> 3           w       NA        NA       NA        NA
```

[`bpe_explore_subsets()`](https://koren6684.github.io/spliv/reference/bpe_explore_subsets.md)
is available for transparent exploratory diagnostics, but exploratory
subgroup search and selecting the first passing rule are not
confirmatory BPE. Confirmatory claims require a pre-specified design and
a reported validation record.
