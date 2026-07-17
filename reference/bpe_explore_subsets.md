# Explore Candidate BPE Subsets

Explores user-supplied candidate subset definitions and reports
diagnostics for transparent exploratory work. This function is not valid
confirmatory BPE by itself. Confirmatory BPE requires a pre-specified
[`bpe_design()`](https://koren6684.github.io/spliv/reference/bpe_design.md)
object supplied to
[`spliv()`](https://koren6684.github.io/spliv/reference/spliv.md) or
[`bpe_validate_design()`](https://koren6684.github.io/spliv/reference/bpe_validate_design.md).

## Usage

``` r
bpe_explore_subsets(
  data,
  spec,
  rules = NULL,
  seed = 1,
  min_n_S = NULL,
  max_F_S = NULL,
  min_varZ_S = NULL,
  return_all = TRUE
)
```

## Arguments

- data:

  Data frame used for estimation.

- spec:

  Either an IV formula (`y ~ X | Z`) or a list with at least `formula`,
  and optional `fe`, `fe_engine`, and `z_names`.

- rules:

  Candidate subset definitions. Each candidate may be a one-sided
  formula evaluated in `data`, a `function(data)` returning a logical
  vector, a character string naming a logical column in `data`, a
  logical vector of length `nrow(data)`, or a named list with components
  `name` and `subset`.

- seed:

  Integer seed for deterministic rule evaluation.

- min_n_S:

  Optional exploratory screen for subset size.

- max_F_S:

  Deprecated exploratory screen for the first-stage F-statistic. The
  first-stage F-statistic is reported as a diagnostic only.

- min_varZ_S:

  Optional exploratory screen for within-subset residualized instrument
  variance.

- return_all:

  If `TRUE` (default), returns all candidates and diagnostics. If
  `FALSE`, returns the first candidate in the user-supplied order. No
  automatic subset search or confirmatory selection is performed.

## Value

If `return_all = TRUE`, a data frame with one row per rule and a list
column `subset`. If `return_all = FALSE`, a list with the first supplied
`subset`, `rule`, and `diagnostics`.

## Details

Warning: this function is exploratory. It should not be used for
confirmatory BPE inference. Confirmatory BPE requires a pre-specified
[`bpe_design()`](https://koren6684.github.io/spliv/reference/bpe_design.md)
object with a substantive rationale.

## Examples

``` r
d <- data.frame(y = rnorm(40), x = rnorm(40), z = rnorm(40),
  inactive = rep(c(TRUE, FALSE), each = 20))
suppressWarnings(bpe_explore_subsets(
  d, y ~ x | z, rules = list(inactive = ~ inactive)))
#>       rule candidate_type available n_S share_S   varZ_S       F_S screen_n_ok
#> 1 inactive        formula      TRUE  20     0.5 1.027887 0.9819404          NA
#>   screen_varZ_ok screen_F_ok message       subset
#> 1             NA          NA         TRUE, TR....
```
