# Extract Tipping-Point Delta from a Sensitivity Path

Returns the tipping-point attribute stored on a
[`spliv_sensitivity_path()`](https://koren6684.github.io/spliv/reference/spliv_sensitivity_path.md)
object. For each term, the tipping point is the smallest `delta` on the
supplied grid whose interval includes zero.

## Usage

``` r
spliv_tipping_point(x)
```

## Arguments

- x:

  A `spliv_sensitivity_path` object.

## Value

Named numeric vector of tipping-point values.

## Examples

``` r
set.seed(6)
d <- data.frame(y = rnorm(80), x = rnorm(80), z = rnorm(80))
p <- spliv_sensitivity_path(y ~ x | z, d, method = "uci", delta_grid = c(0, 0.1))
spliv_tipping_point(p)
#> (Intercept)           x 
#>           0           0 
#> attr(,"message")
#>                                (Intercept) 
#> "Baseline interval already includes zero." 
#>                                          x 
#> "Baseline interval already includes zero." 
```
