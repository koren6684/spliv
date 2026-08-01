# API cleanup audit (0.2.0)

Version 0.2.0 exposes the canonical `spliv*()`/`bpe_*()` workflow together
with genuinely distinct lower-level `sp_*()` sensitivity functions and
`plot_sp_sensitivity()`. Accidental aliases present in the 0.1.x interface
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

`.demean_fixest` and `.demean_lfe` remain internal FE engines. The old prior,
subset, matrix-alias, and parser helpers were removed because they have no
package-internal callers; the named prior embedding helper used by BPE remains
internal.

The requested dead-code review produced these decisions:

- `conley_ltz_mats` and `conley_uci_mats` were not preserved as historical
  names. Their necessary numerical engines remain as the internal
  `.sp_ltz_mats` and `.sp_uci_mats` functions.
- `.resolve_bpe_omega` was removed with the obsolete heuristic/none BPE
  covariance route.
- `.estimate_gamma_zero_first_stage` was removed because no maintained caller
  used it; confirmatory BPE uses the current reduced-form diagnostic engine.
- `.embed_prior_into_full_Z` and `bpe_prior_mats` were removed. The smaller
  `.embed_prior_by_names` helper remains because final BPE estimation calls it.
- `.iv_inst_names` was removed because parsing already records instrument
  names and no maintained caller needed the duplicate helper.

## Fit classes

New estimator results have class `spliv_fit`. Unexported S3 methods for the
historical `plausexog_fit` class remain only to support printing and plotting
objects saved by older package versions.
