# LTZ Sensitivity over Delta Grid

LTZ Sensitivity over Delta Grid

## Usage

``` r
sp_sensitivity_ltz_normal(
  formula,
  data = NULL,
  term,
  inst_vary,
  delta_grid,
  mean_fun = function(delta) 0,
  sd_fun = function(delta) delta,
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

- mean_fun:

  Function mapping delta to prior mean.

- sd_fun:

  Function mapping delta to prior sd.

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
set.seed(10)
d <- data.frame(y = rnorm(60), x = rnorm(60), z = rnorm(60))
sp_sensitivity_ltz_normal(y ~ x | z, d, term = "x", inst_vary = "z",
  delta_grid = c(0, 0.1), scale_instrument = "none")
#>   delta estimate  conf.low conf.high             method
#> 1   0.0 3.154325 -10.85246  17.16111 LTZ (Normal prior)
#> 2   0.1 3.154325 -11.21696  17.52561 LTZ (Normal prior)
```
