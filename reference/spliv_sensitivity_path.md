# Sensitivity Path over Delta Grid

Runs [`spliv()`](https://koren6684.github.io/spliv/reference/spliv.md)
repeatedly over a user-supplied `delta_grid` and returns a tidy
sensitivity-path object. This is the recommended workflow for patterned
or uniform sensitivity analysis: users should usually report how
intervals change over a range of `delta` values rather than selecting
one arbitrary sensitivity level.

## Usage

``` r
spliv_sensitivity_path(
  formula,
  data,
  method = c("uci", "ltz"),
  delta_grid = seq(0, 0.2, by = 0.01),
  violation_pattern = NULL,
  stop_on_error = TRUE,
  ...
)
```

## Arguments

- formula:

  IV formula `y ~ X | Z`.

- data:

  Data frame.

- method:

  One of `"uci"` or `"ltz"`. `spliv_sensitivity_path()` intentionally
  excludes confirmatory BPE.

- delta_grid:

  Non-negative sensitivity grid. For UCI, each value `d` implies theta
  bounds `[-d, +d]`.

- violation_pattern:

  Optional
  [`spliv_pattern()`](https://koren6684.github.io/spliv/reference/spliv_pattern.md)
  object. If omitted, the path uses the package's backward-compatible
  uniform direct-effect pattern.

- stop_on_error:

  Logical; if `TRUE` (default), stop on the first failed fit. If
  `FALSE`, record `NA` rows and store the error message in the returned
  `error` column.

- ...:

  Additional named arguments passed through to
  [`spliv()`](https://koren6684.github.io/spliv/reference/spliv.md),
  such as `fe`, `vcov`, `cluster`, `scale_instrument`, or
  `grid = list(level = 0.95)`.

## Value

A data frame with class `c("spliv_sensitivity_path", "data.frame")`. The
returned object includes path columns such as `delta`, `method`,
`estimate`, `conf_low`, `conf_high`, `contains_zero`, and pattern
metadata, plus attributes containing the original call, the supplied
grid, and a tipping-point summary.

## Examples

``` r
set.seed(5)
d <- data.frame(y = rnorm(80), x = rnorm(80), z = rnorm(80))
p <- spliv_sensitivity_path(y ~ x | z, d, method = "uci",
  delta_grid = c(0, 0.1), vcov = "hc1")
head(p)
#>          term delta method  estimate    conf_low conf_high contains_zero
#> 1 (Intercept)   0.0    uci  0.108665  -0.4663013 0.6836313          TRUE
#> 2           x   0.0    uci -1.979595  -6.9028208 2.9436307          TRUE
#> 3 (Intercept)   0.1    uci  0.108665  -0.6973215 1.0299906          TRUE
#> 4           x   0.1    uci -1.979595 -10.4825196 4.4046121          TRUE
#>            pattern_name pattern_type violation_pattern_used scale_instrument
#> 1 Uniform direct effect      uniform                  FALSE      residual_sd
#> 2 Uniform direct effect      uniform                  FALSE      residual_sd
#> 3 Uniform direct effect      uniform                  FALSE      residual_sd
#> 4 Uniform direct effect      uniform                  FALSE      residual_sd
#>   nobs se theta_min theta_max baseline_estimate baseline_conf_low
#> 1   80 NA       0.0       0.0          0.108665        -0.4663013
#> 2   80 NA       0.0       0.0         -1.979595        -6.9028208
#> 3   80 NA      -0.1       0.1          0.108665        -0.4663013
#> 4   80 NA      -0.1       0.1         -1.979595        -6.9028208
#>   baseline_conf_high crosses_baseline_sign significant_at_level error
#> 1          0.6836313                  TRUE                FALSE  <NA>
#> 2          2.9436307                  TRUE                FALSE  <NA>
#> 3          0.6836313                  TRUE                FALSE  <NA>
#> 4          2.9436307                  TRUE                FALSE  <NA>
```
