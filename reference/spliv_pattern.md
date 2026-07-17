# Create a Patterned Exclusion-Violation Object

Creates a researcher-specified direct-effect pattern for patterned
sensitivity analysis. The pattern determines where direct effects of the
instrument are expected to be larger or smaller in the estimation
sample.

## Usage

``` r
spliv_pattern(
  name,
  pattern,
  rationale,
  variables_used = NULL,
  pattern_type = NULL,
  normalize = c("max_abs", "sd", "none"),
  center = FALSE,
  pre_specified = TRUE,
  notes = NULL
)
```

## Arguments

- name:

  Short pattern name.

- pattern:

  Pattern specification. Accepts a one-sided formula, a
  `function(data)`, a numeric vector of length `nrow(data)`, a logical
  vector of length `nrow(data)`, or a character string naming a numeric
  or logical column in `data`.

- rationale:

  Substantive justification for the proposed pattern.

- variables_used:

  Optional character vector describing the variables used to define the
  pattern.

- pattern_type:

  Optional short label such as `"theory_defined"` or `"design_based"`.

- normalize:

  One of `"max_abs"` (default), `"sd"`, or `"none"`.

- center:

  Logical; if `TRUE`, subtracts the estimation-sample mean before
  normalization.

- pre_specified:

  Logical indicator for whether the pattern was chosen before examining
  outcome results.

- notes:

  Optional free-form notes.

## Value

An object of class `"spliv_pattern"`.

## Examples

``` r
p <- spliv_pattern(
  name = "Exposure", pattern = ~ exposure,
  rationale = "The alternative channel follows exposure.",
  variables_used = "exposure", pattern_type = "theory_defined"
)
p$name
#> [1] "Exposure"
```
