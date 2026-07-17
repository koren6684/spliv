# Demean with fixest

Demeans vectors/matrices by high-dimensional fixed effects using
[`fixest::demean`](https://lrberge.github.io/fixest/reference/demean.html).

## Usage

``` r
demean_fixest(y, X, Z, W = NULL, fe_fml, data)
```

## Arguments

- y:

  Outcome vector.

- X:

  Regressor matrix.

- Z:

  Instrument matrix.

- W:

  Optional controls matrix.

- fe_fml:

  One-sided FE formula.

- data:

  Data frame aligned with the rows of `y`, `X`, `Z`, `W`.

## Value

List with demeaned `y`, `X`, `Z`, `W`, and FE frame.

## Examples

``` r
d <- data.frame(g = factor(rep(1:2, each = 3)))
y <- 1:6; X <- matrix(rnorm(6), 6, 1); Z <- matrix(rnorm(6), 6, 1)
demean_fixest(y, X, Z, fe_fml = ~ g, data = d)
#> $y
#> [1] -1  0  1 -1  0  1
#> 
#> $X
#>            [,1]
#> [1,] -1.4082679
#> [2,]  1.5619304
#> [3,] -0.1536625
#> [4,]  0.2951044
#> [5,]  0.5245657
#> [6,] -0.8196701
#> 
#> $Z
#>            [,1]
#> [1,]  0.1329394
#> [2,]  0.3115066
#> [3,] -0.4444460
#> [4,] -0.3374741
#> [5,]  0.6220733
#> [6,] -0.2845992
#> 
#> $W
#> NULL
#> 
#> $fe_df
#>   g
#> 1 1
#> 2 1
#> 3 1
#> 4 2
#> 5 2
#> 6 2
#> 
```
