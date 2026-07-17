# Build LTZ Prior Matrices for Chosen Instruments

Build LTZ Prior Matrices for Chosen Instruments

## Usage

``` r
sp_prior_ltz(
  formula,
  data = NULL,
  inst_vary,
  mean = 0,
  sd = 1,
  vcov = c("hc1", "hc0", "iid", "cluster"),
  cluster = NULL
)

conley_prior_ltz(
  formula,
  data = NULL,
  inst_vary,
  mean = 0,
  sd = 1,
  vcov = c("hc1", "hc0", "iid", "cluster"),
  cluster = NULL
)
```

## Arguments

- formula:

  Formula or `plausexog_fit`.

- data:

  Data frame when `formula` is a formula.

- inst_vary:

  Instrument(s) with non-degenerate prior.

- mean:

  Prior mean(s).

- sd:

  Prior sd(s).

- vcov:

  Vcov type.

- cluster:

  Optional cluster ids.

## Value

List with `mu`, `omega`, and instrument names.

## Examples

``` r
set.seed(9)
d <- data.frame(y = rnorm(60), x = rnorm(60), z = rnorm(60))
sp_prior_ltz(y ~ x | z, d, inst_vary = "z", sd = 0.1)
#> $mu
#> [1] 0 0
#> 
#> $omega
#>      [,1] [,2]
#> [1,]    0 0.00
#> [2,]    0 0.01
#> 
#> $inst_names
#> [1] "(Intercept)" "z"          
#> 
```
