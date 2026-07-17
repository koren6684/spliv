# Create a Confirmatory BPE Design Object

Creates a researcher-specified design object for confirmatory BPE. The
subset must be justified outside the outcome analysis and should
represent a theoretically motivated instrument-inactive subset.

## Usage

``` r
bpe_design(
  name,
  subset,
  rationale,
  variables_used = NULL,
  subset_type = NULL,
  pre_specified = TRUE,
  transportability_rationale = NULL,
  notes = NULL
)
```

## Arguments

- name:

  Short design name.

- subset:

  Subset specification. Accepts a one-sided formula, a `function(data)`,
  a logical vector of length `nrow(data)`, or a character string naming
  a logical column in `data`.

- rationale:

  Non-empty substantive justification for the subset.

- variables_used:

  Optional character vector describing the design variables used to
  define the subset. This is required for function-based subsets unless
  the package can safely infer the variables.

- subset_type:

  Optional short label such as `"theory_defined"` or `"design_based"`.

- pre_specified:

  Logical indicator for confirmatory use. Confirmatory BPE requires
  `TRUE`.

- transportability_rationale:

  Optional description of why the direct effect learned in the subset
  may be informative for the target sample.

- notes:

  Optional free-form notes.

## Value

An object of class `"spliv_bpe_design"`.

## Examples

``` r
design <- bpe_design(
  "Inactive subset", ~ inactive,
  rationale = "The treatment channel is absent in this subset.",
  variables_used = "inactive", pre_specified = TRUE
)
bpe_eval_subset(design, data.frame(inactive = c(TRUE, FALSE)))
#> [1]  TRUE FALSE
#> attr(,"bpe_warnings")
#> character(0)
#> attr(,"bpe_na_share")
#> [1] 0
```
