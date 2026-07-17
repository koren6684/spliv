# Demean with lfe

Demeans vectors/matrices by high-dimensional fixed effects using
[`lfe::demeanlist`](https://rdrr.io/pkg/lfe/man/demeanlist.html).

## Usage

``` r
demean_lfe(y, X, Z, W = NULL, fe_list)
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

- fe_list:

  List/data.frame of FE ids.

## Value

List with demeaned `y`, `X`, `Z`, `W`.

## Examples

``` r
d <- data.frame(g = factor(rep(1:2, each = 3)))
y <- 1:6; X <- matrix(rnorm(6), 6, 1); Z <- matrix(rnorm(6), 6, 1)
demean_lfe(y, X, Z, fe_list = d["g"])
#> $y
#> [1] -1  0  1 -1  0  1
#> 
#> $X
#>            [,1]
#> [1,] -0.2928923
#> [2,] -0.4350513
#> [3,]  0.7279436
#> [4,]  0.5186932
#> [5,] -0.2686336
#> [6,] -0.2500596
#> 
#> $Z
#>            [,1]
#> [1,]  0.2630599
#> [2,] -0.8886963
#> [3,]  0.6256364
#> [4,]  0.6800044
#> [5,]  0.4563673
#> [6,] -1.1363717
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
