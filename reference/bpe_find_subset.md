# Deprecated Exploratory BPE Subset Search

`bpe_find_subset()` is retained for backward compatibility only. Use
[`bpe_explore_subsets()`](https://koren6684.github.io/spliv/reference/bpe_explore_subsets.md)
for exploratory work and convert any theory-justified subset into a
confirmatory
[`bpe_design()`](https://koren6684.github.io/spliv/reference/bpe_design.md)
object before calling
[`spliv()`](https://koren6684.github.io/spliv/reference/spliv.md).

## Usage

``` r
bpe_find_subset(
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

See
[`bpe_explore_subsets()`](https://koren6684.github.io/spliv/reference/bpe_explore_subsets.md).

## Examples

``` r
d <- data.frame(y = rnorm(40), x = rnorm(40), z = rnorm(40),
  inactive = rep(c(TRUE, FALSE), each = 20))
suppressWarnings(bpe_find_subset(
  d, y ~ x | z, rules = list(inactive = ~ inactive)))
#>       rule candidate_type available n_S share_S    varZ_S       F_S screen_n_ok
#> 1 inactive        formula      TRUE  20     0.5 0.8294879 0.3497427          NA
#>   screen_varZ_ok screen_F_ok message       subset
#> 1             NA          NA         TRUE, TR....
```
