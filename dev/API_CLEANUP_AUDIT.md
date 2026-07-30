# API cleanup audit (0.2.0)

Version 0.2.0 exposes the canonical `spliv*()`/`bpe_*()` workflow together
with genuinely distinct lower-level `sp_*()` sensitivity functions and
`plot_sp_sensitivity()`. Accidental aliases from the 0.1.1 staging history
were removed rather than retained as a deprecation transition.

## Canonical public API

`spliv`, `spliv_pattern`, `spliv_eval_pattern`, `spliv_sensitivity_path`,
`spliv_tipping_point`, `bpe_design`, `bpe_eval_subset`,
`bpe_validate_design`, and `bpe_explore_subsets` are the public workflow.
Exploratory subset search remains explicitly non-confirmatory BPE.

## Advanced low-level interface

The exported low-level functions are `sp_ltz`, `sp_uci`, `sp_prior_ltz`,
`sp_sensitivity_ltz_normal`, `sp_sensitivity_uci_support`,
`sp_sensitivity_ltz_uniform01_as_normal`, and `plot_sp_sensitivity`.
They expose distinct prior, bound, and plotting controls; ordinary users
should normally start with `spliv()` and `spliv_sensitivity_path()`.

## Internal helpers

`.demean_fixest`, `.demean_lfe`, `.embed_prior_into_full_Z`, and
`.estimate_gamma_zero_first_stage` remain unexported implementation helpers.
`iv_inst_names` was also made internal because it has no package-internal
callers.

## Fit classes

New estimator results have class `spliv_fit`. Unexported S3 methods for the
historical `plausexog_fit` class remain only to support printing and plotting
objects saved by older package versions.
