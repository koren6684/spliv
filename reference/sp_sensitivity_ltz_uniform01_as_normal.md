# LTZ Sensitivity with Normal Approximation to U(0, delta)

LTZ Sensitivity with Normal Approximation to U(0, delta)

## Usage

``` r
sp_sensitivity_ltz_uniform01_as_normal(
  formula,
  data = NULL,
  term,
  inst_vary,
  delta_grid,
  level = 0.95,
  vcov = c("hc1", "hc0", "iid", "cluster"),
  cluster = NULL,
  scale_instrument = c("residual_sd", "none")
)
```

## Arguments

- formula:

  Formula or `spliv_fit`.

- data:

  Data frame when `formula` is a formula.

- term:

  Coefficient name to track.

- inst_vary:

  Instrument(s) with plausible violation.

- delta_grid:

  Delta grid.

- level:

  Confidence level.

- vcov:

  Vcov type.

- cluster:

  Optional cluster ids.

- scale_instrument:

  One of `"residual_sd"` (default) or `"none"`.

## Value

Data frame with sensitivity path.

## Examples

``` r
set.seed(12)
d <- data.frame(y = rnorm(60), x = rnorm(60), z = rnorm(60))
sp_sensitivity_ltz_uniform01_as_normal(y ~ x | z, d, term = "x",
  inst_vary = "z", delta_grid = c(0, 0.1), scale_instrument = "none")
#>   delta    estimate   conf.low conf.high                            method
#> 1   0.0 -0.03735536 -1.1100982  1.035387 LTZ (Normal approx to U(0,delta))
#> 2   0.1  0.26564199 -0.8605619  1.391846 LTZ (Normal approx to U(0,delta))
```
