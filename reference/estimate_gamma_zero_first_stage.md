# Legacy Exploratory BPE Prior Helper

`estimate_gamma_zero_first_stage()` is retained as a legacy exploratory
helper. It is not sufficient for confirmatory BPE inference.
Confirmatory BPE requires
[`bpe_design()`](https://koren6684.github.io/spliv/reference/bpe_design.md),
[`bpe_validate_design()`](https://koren6684.github.io/spliv/reference/bpe_validate_design.md),
and `spliv(method = "bpe", ...)`.

## Usage

``` r
estimate_gamma_zero_first_stage(
  data,
  y_name,
  z_names,
  controls = character(0),
  subset,
  fe = NULL
)
```

## Arguments

- data:

  Data frame.

- y_name:

  Outcome variable name.

- z_names:

  Instrument names whose direct effects are estimated.

- controls:

  Optional character vector of additional controls.

- subset:

  Legacy subset specification. Accepts a
  [`bpe_design()`](https://koren6684.github.io/spliv/reference/bpe_design.md)
  object, a one-sided formula, a `function(data)`, a logical vector, or
  a character string naming a logical column in `data`.

- fe:

  Optional one-sided FE formula. If supplied, uses
  [`fixest::feols`](https://lrberge.github.io/fixest/reference/feols.html).

## Value

List with `mu_hat`, `omega_hat`, and diagnostics.

## Examples

``` r
d <- data.frame(y = rnorm(40), z = rnorm(40), inactive = rep(c(TRUE, FALSE), each = 20))
suppressWarnings(estimate_gamma_zero_first_stage(
  d, y_name = "y", z_names = "z", subset = ~ inactive))
#> $mu_hat
#> [1] 0.0788356
#> 
#> $omega_hat
#>        [,1]
#> z 0.0786822
#> 
#> $diagnostics
#> $diagnostics$n_subset
#> [1] 20
#> 
#> $diagnostics$p
#> [1] 2
#> 
#> $diagnostics$engine
#> [1] "ols"
#> 
#> 
```
