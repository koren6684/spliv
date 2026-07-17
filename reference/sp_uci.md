# Union of Confidence Intervals for Plausibly Exogenous IV

Formula-level wrapper over the matrix core UCI implementation.

## Usage

``` r
sp_uci(
  formula,
  data,
  inst,
  gmin,
  gmax,
  grid = 21,
  level = 0.95,
  vcov = c("hc1", "hc0", "iid", "cluster"),
  cluster = NULL
)

conley_uci(
  formula,
  data,
  inst,
  gmin,
  gmax,
  grid = 21,
  level = 0.95,
  vcov = c("hc1", "hc0", "iid", "cluster"),
  cluster = NULL
)
```

## Arguments

- formula:

  IV formula `y ~ X | Z`.

- data:

  Data frame.

- inst:

  Instrument names to vary.

- gmin:

  Lower bound(s) on gamma.

- gmax:

  Upper bound(s) on gamma.

- grid:

  Grid size per varied instrument.

- level:

  Confidence level.

- vcov:

  One of `"iid"`, `"hc0"`, `"hc1"`, `"cluster"`.

- cluster:

  Cluster ids when `vcov = "cluster"`.

## Value

Data frame with union confidence intervals.

## Examples

``` r
set.seed(4)
d <- data.frame(y = rnorm(80), x = rnorm(80), z = rnorm(80))
sp_uci(y ~ x | z, d, inst = "z", gmin = -0.1, gmax = 0.1, grid = 5)
#>          term   conf.low conf.high
#> 1 (Intercept) -0.4302267 0.5464528
#> 2           x -5.4262559 3.0779421
```
