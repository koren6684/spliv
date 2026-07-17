# Local-to-Zero Inference for Plausibly Exogenous IV

Formula-level wrapper over the matrix core LTZ implementation.

## Usage

``` r
sp_ltz(
  formula,
  data,
  omega,
  mu,
  level = 0.95,
  vcov = c("hc1", "hc0", "iid", "cluster"),
  cluster = NULL
)

conley_ltz(
  formula,
  data,
  omega,
  mu,
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

- omega:

  Prior covariance matrix over instrument direct effects.

- mu:

  Prior mean vector over instrument direct effects.

- level:

  Confidence level.

- vcov:

  One of `"iid"`, `"hc0"`, `"hc1"`, `"cluster"`.

- cluster:

  Cluster ids when `vcov = "cluster"`.

## Value

Data frame with adjusted estimates and confidence intervals.

## Examples

``` r
set.seed(3)
d <- data.frame(y = rnorm(80), x = rnorm(80), z = rnorm(80))
prior <- conley_prior_ltz(y ~ x | z, d, inst_vary = "z", sd = 0.1)
sp_ltz(y ~ x | z, d, prior$omega, prior$mu)
#>          term     estimate std.error   conf.low  conf.high
#> 1 (Intercept) -0.003841335  0.175378 -0.3475759  0.3398932
#> 2           x  1.221902261  5.489363 -9.5370506 11.9808551
```
