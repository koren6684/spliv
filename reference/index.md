# Package index

## Get started

- [`spliv-package`](https://koren6684.github.io/spliv/reference/spliv-package.md)
  : spliv: Patterned Sensitivity Analysis for IV

## Main estimation

- [`plausexog()`](https://koren6684.github.io/spliv/reference/spliv.md)
  [`plausexog_iv()`](https://koren6684.github.io/spliv/reference/spliv.md)
  [`spliv()`](https://koren6684.github.io/spliv/reference/spliv.md) :
  Patterned Sensitivity Analysis for Plausibly Exogenous IV
- [`sp_ltz()`](https://koren6684.github.io/spliv/reference/sp_ltz.md)
  [`conley_ltz()`](https://koren6684.github.io/spliv/reference/sp_ltz.md)
  : Local-to-Zero Inference for Plausibly Exogenous IV
- [`sp_uci()`](https://koren6684.github.io/spliv/reference/sp_uci.md)
  [`conley_uci()`](https://koren6684.github.io/spliv/reference/sp_uci.md)
  : Union of Confidence Intervals for Plausibly Exogenous IV

## Pattern definitions

- [`spliv_pattern()`](https://koren6684.github.io/spliv/reference/spliv_pattern.md)
  : Create a Patterned Exclusion-Violation Object
- [`spliv_eval_pattern()`](https://koren6684.github.io/spliv/reference/spliv_eval_pattern.md)
  : Evaluate a Direct-Effect Pattern
- [`iv_inst_names()`](https://koren6684.github.io/spliv/reference/iv_inst_names.md)
  : Instrument Names for Parsed IV Formula

## Sensitivity paths

- [`spliv_sensitivity_path()`](https://koren6684.github.io/spliv/reference/spliv_sensitivity_path.md)
  : Sensitivity Path over Delta Grid
- [`spliv_tipping_point()`](https://koren6684.github.io/spliv/reference/spliv_tipping_point.md)
  : Extract Tipping-Point Delta from a Sensitivity Path
- [`sp_sensitivity_ltz_normal()`](https://koren6684.github.io/spliv/reference/sp_sensitivity_ltz_normal.md)
  [`conley_sensitivity_ltz_normal()`](https://koren6684.github.io/spliv/reference/sp_sensitivity_ltz_normal.md)
  : LTZ Sensitivity over Delta Grid
- [`sp_sensitivity_uci_support()`](https://koren6684.github.io/spliv/reference/sp_sensitivity_uci_support.md)
  [`conley_sensitivity_uci_support()`](https://koren6684.github.io/spliv/reference/sp_sensitivity_uci_support.md)
  : UCI Sensitivity over Delta Grid

## BPE design and validation

- [`bpe_design()`](https://koren6684.github.io/spliv/reference/bpe_design.md)
  : Create a Confirmatory BPE Design Object
- [`bpe_eval_subset()`](https://koren6684.github.io/spliv/reference/bpe_eval_subset.md)
  : Evaluate a Confirmatory BPE Subset
- [`bpe_validate_design()`](https://koren6684.github.io/spliv/reference/bpe_validate_design.md)
  : Validate a Confirmatory BPE Design
- [`embed_prior_into_full_Z()`](https://koren6684.github.io/spliv/reference/embed_prior_into_full_Z.md)
  : Embed Prior into Full Instrument Space

## Plotting

- [`plot_sp_sensitivity()`](https://koren6684.github.io/spliv/reference/plot_sp_sensitivity.md)
  [`plot_conley_sensitivity()`](https://koren6684.github.io/spliv/reference/plot_sp_sensitivity.md)
  : Plot Patterned Sensitivity Output or a Fitted Object

## Compatibility/deprecated functions

- [`plausexog()`](https://koren6684.github.io/spliv/reference/spliv.md)
  [`plausexog_iv()`](https://koren6684.github.io/spliv/reference/spliv.md)
  [`spliv()`](https://koren6684.github.io/spliv/reference/spliv.md) :
  Patterned Sensitivity Analysis for Plausibly Exogenous IV
- [`bpe_explore_subsets()`](https://koren6684.github.io/spliv/reference/bpe_explore_subsets.md)
  : Explore Candidate BPE Subsets
- [`bpe_find_subset()`](https://koren6684.github.io/spliv/reference/bpe_find_subset.md)
  : Deprecated Exploratory BPE Subset Search
- [`estimate_gamma_zero_first_stage()`](https://koren6684.github.io/spliv/reference/estimate_gamma_zero_first_stage.md)
  : Legacy Exploratory BPE Prior Helper

## Compatibility and legacy helpers

Compatibility functions and lower-level helpers retained for existing
workflows.

- [`demean_fixest()`](https://koren6684.github.io/spliv/reference/demean_fixest.md)
  : Demean with fixest
- [`demean_lfe()`](https://koren6684.github.io/spliv/reference/demean_lfe.md)
  : Demean with lfe
- [`sp_prior_ltz()`](https://koren6684.github.io/spliv/reference/sp_prior_ltz.md)
  [`conley_prior_ltz()`](https://koren6684.github.io/spliv/reference/sp_prior_ltz.md)
  : Build LTZ Prior Matrices for Chosen Instruments
- [`sp_sensitivity_ltz_uniform01_as_normal()`](https://koren6684.github.io/spliv/reference/sp_sensitivity_ltz_uniform01_as_normal.md)
  [`conley_sensitivity_ltz_uniform01_as_normal()`](https://koren6684.github.io/spliv/reference/sp_sensitivity_ltz_uniform01_as_normal.md)
  : LTZ Sensitivity with Normal Approximation to U(0, delta)
