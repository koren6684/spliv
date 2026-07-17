# UCI Sensitivity over Delta Grid

UCI Sensitivity over Delta Grid

## Usage

``` r
sp_sensitivity_uci_support(
  formula,
  data = NULL,
  term,
  inst_vary,
  delta_grid,
  gmin_fun = function(delta) -delta,
  gmax_fun = function(delta) delta,
  grid = 41,
  level = 0.95,
  vcov = c("hc1", "hc0", "iid", "cluster"),
  cluster = NULL,
  scale_instrument = c("residual_sd", "none")
)

conley_sensitivity_uci_support(
  formula,
  data = NULL,
  term,
  inst_vary,
  delta_grid,
  gmin_fun = function(delta) -delta,
  gmax_fun = function(delta) delta,
  grid = 41,
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

  Instrument(s) whose gamma bounds vary.

- delta_grid:

  Delta grid.

- gmin_fun:

  Function mapping delta to lower bound.

- gmax_fun:

  Function mapping delta to upper bound.

- grid:

  Points per instrument bound.

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
