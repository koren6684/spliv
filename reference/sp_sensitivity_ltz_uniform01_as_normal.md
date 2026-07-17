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

conley_sensitivity_ltz_uniform01_as_normal(
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

  Formula or `plausexog_fit`.

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
