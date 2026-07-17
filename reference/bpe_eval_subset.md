# Evaluate a Confirmatory BPE Subset

Evaluates a
[`bpe_design()`](https://koren6684.github.io/spliv/reference/bpe_design.md)
object on a data frame and returns the resulting logical indicator.

## Usage

``` r
bpe_eval_subset(design, data, max_na_share = 0.05)
```

## Arguments

- design:

  A
  [`bpe_design()`](https://koren6684.github.io/spliv/reference/bpe_design.md)
  object.

- data:

  Data frame used for evaluation.

- max_na_share:

  Maximum allowable share of `NA` values in the subset indicator before
  evaluation fails. The default is `0.05`.

## Value

A logical vector of length `nrow(data)`.

## Examples

``` r
design <- bpe_design("Inactive", ~ inactive,
  rationale = "The treatment channel is absent.")
bpe_eval_subset(design, data.frame(inactive = c(TRUE, FALSE)))
#> [1]  TRUE FALSE
#> attr(,"bpe_warnings")
#> character(0)
#> attr(,"bpe_na_share")
#> [1] 0
```
