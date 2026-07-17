# Getting started with spliv

## Baseline IV and uniform sensitivity

The [`spliv()`](https://koren6684.github.io/spliv/reference/spliv.md)
estimator accepts the usual IV formula `y ~ X | Z`. A baseline fit with
`delta = 0` is an ordinary IV calculation, while a positive `delta`
allows a bounded uniform direct effect.

``` r

baseline <- spliv(f, d, method = "uci", delta = 0, vcov = "hc1")
baseline$estimates
#>          term    conf.low  conf.high
#> 1 (Intercept) -0.18384632 0.09088447
#> 2           x  1.29870349 1.84425108
#> 3           w -0.01990612 0.28613695
uniform <- spliv(f, d, method = "uci", delta = 0.20, vcov = "hc1",
                 grid = list(steps = 9))
uniform$estimates
#>          term   conf.low conf.high
#> 1 (Intercept) -0.2487226 0.1183848
#> 2           x  0.8961868 2.3416858
#> 3           w -0.2070957 0.4255989
```

## A patterned UCI/LTZ analysis

Patterns are explicit objects with a substantive rationale. Here the
direct effect is allowed to be larger where the synthetic exposure `w`
is larger.

``` r

pattern <- spliv_pattern(
  name = "Exposure pattern", pattern = ~ w,
  rationale = "The alternative channel is stronger at higher exposure.",
  variables_used = "w", pattern_type = "theory_defined",
  normalize = "max_abs"
)
spliv_eval_pattern(pattern, d)[1:5]
#> [1]  0.083746179  0.035128314  0.380408513  0.233287898 -0.001435305
patterned_uci <- spliv(f, d, method = "uci", delta = 0.20, vcov = "hc1",
                       violation_pattern = pattern, grid = list(steps = 9))
patterned_ltz <- spliv(f, d, method = "ltz", delta = 0.20, vcov = "hc1",
                       violation_pattern = pattern)
patterned_uci$estimates
#>          term    conf.low  conf.high
#> 1 (Intercept) -0.18444170 0.09195752
#> 2           x  1.29090596 1.85535775
#> 3           w -0.02698832 0.30162026
patterned_ltz$estimates
#>          term    estimate  std.error    conf.low  conf.high
#> 1 (Intercept) -0.04648093 0.07008608 -0.18384711 0.09088526
#> 2           x  1.57147728 0.13918270  1.29868421 1.84427036
#> 3           w  0.13311542 0.07874347 -0.02121894 0.28744978
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
#>          term delta method    estimate    conf_low  conf_high contains_zero
#> 1 (Intercept)  0.00    uci -0.04648093 -0.18384632 0.09088447          TRUE
#> 2           x  0.00    uci  1.57147728  1.29870349 1.84425108         FALSE
#> 3           w  0.00    uci  0.13311542 -0.01990612 0.28613695          TRUE
#> 4 (Intercept)  0.05    uci -0.04648093 -0.18393310 0.09109067          TRUE
#> 5           x  0.05    uci  1.57147728  1.29713579 1.84664606         FALSE
#> 6           w  0.05    uci  0.13311542 -0.02147632 0.28982073          TRUE
#>       pattern_name   pattern_type violation_pattern_used scale_instrument nobs
#> 1 Exposure pattern theory_defined                   TRUE      residual_sd  240
#> 2 Exposure pattern theory_defined                   TRUE      residual_sd  240
#> 3 Exposure pattern theory_defined                   TRUE      residual_sd  240
#> 4 Exposure pattern theory_defined                   TRUE      residual_sd  240
#> 5 Exposure pattern theory_defined                   TRUE      residual_sd  240
#> 6 Exposure pattern theory_defined                   TRUE      residual_sd  240
#>   se theta_min theta_max baseline_estimate baseline_conf_low baseline_conf_high
#> 1 NA      0.00      0.00       -0.04648093       -0.18384632         0.09088447
#> 2 NA      0.00      0.00        1.57147728        1.29870349         1.84425108
#> 3 NA      0.00      0.00        0.13311542       -0.01990612         0.28613695
#> 4 NA     -0.05      0.05       -0.04648093       -0.18384632         0.09088447
#> 5 NA     -0.05      0.05        1.57147728        1.29870349         1.84425108
#> 6 NA     -0.05      0.05        0.13311542       -0.01990612         0.28613695
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
outcome analysis. Validation reports eligibility diagnostics; it does
not search across candidate subgroups.

``` r

design <- bpe_design(
  name = "Theory-defined inactive subset", subset = ~ inactive,
  rationale = "The treatment channel is absent in the inactive subset.",
  variables_used = "inactive", subset_type = "theory_defined",
  pre_specified = TRUE,
  transportability_rationale = "The subset direct effect is informative for the target sample."
)
validation <- bpe_validate_design(
  f, d, design = design, vcov = "hc1",
  bpe_min_n_S = 40,
  bpe_equiv_margin = 0.25 * sd(resid(lm(x ~ w)))
)
validation[c("n_S", "equivalence_passed", "eligibility_passed")]
#> $n_S
#> [1] 120
#> 
#> $equivalence_passed
#> [1] TRUE
#> 
#> $eligibility_passed
#> [1] TRUE

# Illustrative synthetic scale only. In a substantive analysis, pre-specify
# the equivalence margin from the scientific design rather than tuning it to
# obtain BPE eligibility.
bpe_margin <- 0.25 * sd(resid(lm(x ~ w)))
bpe_fit <- spliv(f, d, method = "bpe", bpe_design = design,
                 vcov = "hc1", bpe_min_n_S = 40,
                 bpe_equiv_margin = bpe_margin)
bpe_fit$estimates
#>          term    estimate  std.error    conf.low conf.high
#> 1 (Intercept) -0.04159845 0.07364906 -0.18594795 0.1027510
#> 2           x  1.51206522 0.30855850  0.90730168 2.1168288
#> 3           w  0.15414021 0.12487183 -0.09060407 0.3988845
```

[`bpe_explore_subsets()`](https://koren6684.github.io/spliv/reference/bpe_explore_subsets.md)
is available for transparent exploratory diagnostics, but exploratory
subgroup search and selecting the first passing rule are not
confirmatory BPE. Confirmatory claims require a pre-specified design and
a reported validation record.
