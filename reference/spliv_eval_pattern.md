# Evaluate a Direct-Effect Pattern

Evaluates a
[`spliv_pattern()`](https://koren6684.github.io/spliv/reference/spliv_pattern.md)
object on a data frame and returns a numeric vector after optional
centering and normalization.

## Usage

``` r
spliv_eval_pattern(pattern, data)
```

## Arguments

- pattern:

  A
  [`spliv_pattern()`](https://koren6684.github.io/spliv/reference/spliv_pattern.md)
  object.

- data:

  Data frame used for evaluation.

## Value

Numeric vector of length `nrow(data)`.

## Examples

``` r
d <- data.frame(exposure = seq(-1, 1, length.out = 5))
p <- spliv_pattern("Exposure", ~ exposure,
  rationale = "Illustration of a monotone exposure pattern.")
spliv_eval_pattern(p, d)
#> [1] -1.0 -0.5  0.0  0.5  1.0
#> attr(,"spliv_pattern_raw_summary")
#> attr(,"spliv_pattern_raw_summary")$n
#> [1] 5
#> 
#> attr(,"spliv_pattern_raw_summary")$mean
#> [1] 0
#> 
#> attr(,"spliv_pattern_raw_summary")$sd
#> [1] 0.7905694
#> 
#> attr(,"spliv_pattern_raw_summary")$min
#> [1] -1
#> 
#> attr(,"spliv_pattern_raw_summary")$max
#> [1] 1
#> 
#> attr(,"spliv_pattern_summary")
#> attr(,"spliv_pattern_summary")$n
#> [1] 5
#> 
#> attr(,"spliv_pattern_summary")$mean
#> [1] 0
#> 
#> attr(,"spliv_pattern_summary")$sd
#> [1] 0.7905694
#> 
#> attr(,"spliv_pattern_summary")$min
#> [1] -1
#> 
#> attr(,"spliv_pattern_summary")$max
#> [1] 1
#> 
#> attr(,"spliv_pattern_warnings")
#> character(0)
```
