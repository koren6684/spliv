# Embed Prior into Full Instrument Space

Expands a prior estimated on a subset of instruments into the full
instrument vector implied by `y ~ X | Z`.

## Usage

``` r
embed_prior_into_full_Z(formula, data, z_names, mu_hat, omega_hat)
```

## Arguments

- formula:

  IV formula `y ~ X | Z`.

- data:

  Data frame.

- z_names:

  Instrument names with prior moments.

- mu_hat:

  Prior means for `z_names`.

- omega_hat:

  Prior covariance for `z_names`.

## Value

List with full-length `mu`, `Omega`, and instrument names.

## Examples

``` r
d <- data.frame(y = rnorm(30), x = rnorm(30), z1 = rnorm(30), z2 = rnorm(30))
embed_prior_into_full_Z(y ~ x | z1 + z2, d,
  z_names = "z1", mu_hat = 0, omega_hat = matrix(0.1, 1, 1))
#> $mu
#> (Intercept)          z1          z2 
#>           0           0           0 
#> 
#> $Omega
#>             (Intercept)  z1 z2
#> (Intercept)           0 0.0  0
#> z1                    0 0.1  0
#> z2                    0 0.0  0
#> 
#> $inst_names
#> [1] "(Intercept)" "z1"          "z2"         
#> 
```
