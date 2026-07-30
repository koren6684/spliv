# Public API audit

This inventory intentionally keeps all existing exports. Categories describe
current positioning, not removal decisions.

## primary_current_api

`spliv`, `sp_ltz`, `sp_uci`, `conley_ltz`, `conley_uci`,
`spliv_pattern`, `spliv_eval_pattern`, `spliv_sensitivity_path`,
`spliv_tipping_point`, `bpe_design`, `bpe_eval_subset`,
`bpe_validate_design`, `plot_sp_sensitivity`, `plot.spliv_sensitivity_path`,
`plot.plausexog_fit`, `plot_conley_sensitivity`, `sp_prior_ltz`,
`sp_sensitivity_ltz_normal`, `sp_sensitivity_uci_support`,
`sp_sensitivity_ltz_uniform01_as_normal`, `iv_inst_names`.

These functions are the documented uniform/patterned sensitivity, path, BPE,
and plotting workflows. Their parameter and return-value documentation is
generated from the roxygen blocks in `R/` and shared aliases are noted in the
generated `.Rd` files.

## compatibility_wrapper

`conley_prior_ltz`, `conley_sensitivity_ltz_normal`,
`conley_sensitivity_uci_support`, `conley_sensitivity_ltz_uniform01_as_normal`,
`plausexog_iv`, `demean_fixest`, `demean_lfe`, `embed_prior_into_full_Z`.

These exports preserve the prior API or provide a thin name-compatible wrapper;
new code should prefer the primary functions described above.

## deprecated_exploratory

`plausexog`, `bpe_explore_subsets`, `bpe_find_subset`,
`estimate_gamma_zero_first_stage`.

They emit explicit deprecation or exploratory warnings where applicable. In
particular, subgroup search is not confirmatory BPE: use a pre-specified
`bpe_design()` and `bpe_validate_design()` for confirmatory work.

## internal_candidate

None of the current exports are removed in this release. Future API review may
consider `demean_*`, `embed_prior_into_full_Z`, or low-level `sp_*` wrappers for
internal status after downstream usage is known; that decision is deferred.

## Documentation status

All exported functions have roxygen parameter/return documentation directly or
through an explicit shared `@rdname` block. Runnable synthetic examples cover
the primary workflows in the package README and the getting-started vignette;
deprecated wrappers are tested for warning behavior. The generated manual is
the source of truth after `roxygen2::roxygenise()`.
