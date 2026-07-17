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

conley_sensitivity_ltz_normal(
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

  Formula or `plausexog_fit`.

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
