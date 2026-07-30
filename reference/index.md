# Package index

## Get started

- [`spliv-package`](https://koren6684.github.io/spliv/reference/spliv-package.md)
  : spliv: Patterned Sensitivity Analysis for IV

## Main estimation

- [`spliv()`](https://koren6684.github.io/spliv/reference/spliv.md) :
  Patterned Sensitivity Analysis for Plausibly Exogenous IV

## Pattern definition

- [`spliv_pattern()`](https://koren6684.github.io/spliv/reference/spliv_pattern.md)
  : Create a Patterned Exclusion-Violation Object
- [`spliv_eval_pattern()`](https://koren6684.github.io/spliv/reference/spliv_eval_pattern.md)
  : Evaluate a Direct-Effect Pattern

## Sensitivity paths and tipping points

- [`spliv_sensitivity_path()`](https://koren6684.github.io/spliv/reference/spliv_sensitivity_path.md)
  : Sensitivity Path over Delta Grid
- [`spliv_tipping_point()`](https://koren6684.github.io/spliv/reference/spliv_tipping_point.md)
  : Extract Tipping-Point Delta from a Sensitivity Path

## Confirmatory BPE

- [`bpe_design()`](https://koren6684.github.io/spliv/reference/bpe_design.md)
  : Create a Confirmatory BPE Design Object
- [`bpe_eval_subset()`](https://koren6684.github.io/spliv/reference/bpe_eval_subset.md)
  : Evaluate a Confirmatory BPE Subset
- [`bpe_validate_design()`](https://koren6684.github.io/spliv/reference/bpe_validate_design.md)
  : Validate a Confirmatory BPE Design

## Exploratory diagnostics

- [`bpe_explore_subsets()`](https://koren6684.github.io/spliv/reference/bpe_explore_subsets.md)
  : Explore Candidate BPE Subsets

## Advanced low-level interface

Lower-level interfaces retained for custom workflows. Ordinary users
should normally prefer
[`spliv()`](https://koren6684.github.io/spliv/reference/spliv.md) and
[`spliv_sensitivity_path()`](https://koren6684.github.io/spliv/reference/spliv_sensitivity_path.md).

- [`sp_ltz()`](https://koren6684.github.io/spliv/reference/sp_ltz.md) :
  Local-to-Zero Inference for Plausibly Exogenous IV
- [`sp_uci()`](https://koren6684.github.io/spliv/reference/sp_uci.md) :
  Union of Confidence Intervals for Plausibly Exogenous IV
- [`sp_prior_ltz()`](https://koren6684.github.io/spliv/reference/sp_prior_ltz.md)
  : Build LTZ Prior Matrices for Chosen Instruments
- [`sp_sensitivity_ltz_normal()`](https://koren6684.github.io/spliv/reference/sp_sensitivity_ltz_normal.md)
  : LTZ Sensitivity over Delta Grid
- [`sp_sensitivity_uci_support()`](https://koren6684.github.io/spliv/reference/sp_sensitivity_uci_support.md)
  : UCI Sensitivity over Delta Grid
- [`sp_sensitivity_ltz_uniform01_as_normal()`](https://koren6684.github.io/spliv/reference/sp_sensitivity_ltz_uniform01_as_normal.md)
  : LTZ Sensitivity with Normal Approximation to U(0, delta)

## Plotting

- [`plot_sp_sensitivity()`](https://koren6684.github.io/spliv/reference/plot_sp_sensitivity.md)
  : Plot Patterned Sensitivity Output or a Fitted Object
