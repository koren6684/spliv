# Baseline before final hardening

Captured 2026-07-31 before the final scientific/API hardening edits. This file
preserves the state against which numerical and reference-output changes are
evaluated.

## Version control and package metadata

- Package branch: `main`
- Package commit: `41993897f83f9388abd2cd13368c2e36ae91aefa`
- Reproducibility branch: `main`
- Reproducibility commit: `27cbb8550210b1f77eb3e3759ba6ff52e47c0303`
- DESCRIPTION version: `0.2.0`
- Reproducibility worktree: clean.
- Package worktree already contained the uncommitted 0.2.0 cleanup changes
  listed below; these are preserved as user-owned work.

Package worktree status before this pass:

```text
 M CITATION.cff
 M NEWS.md
 M R/bpe.R
 M R/core_mats.R
 M R/iv_parse.R
 M R/plotting.R
 M R/spliv.R
 M README.md
 M dev/API_CLEANUP_AUDIT.md
 M man/bpe_explore_subsets.Rd
 M man/plot_sp_sensitivity.Rd
 M man/sp_sensitivity_ltz_normal.Rd
 M man/sp_sensitivity_ltz_uniform01_as_normal.Rd
 M man/sp_sensitivity_uci_support.Rd
 M man/spliv.Rd
 M man/spliv_sensitivity_path.Rd
 M tests/testthat/test-core-api.R
 M vignettes/getting-started.Rmd
?? dev/RJOURNAL_PACKAGE_AUDIT.md
?? dev/benchmarks/
?? tests/testthat/test-rjournal-hardening.R
```

## Public API baseline

Exports (exactly 16):

```text
bpe_design
bpe_eval_subset
bpe_explore_subsets
bpe_validate_design
plot_sp_sensitivity
sp_ltz
sp_prior_ltz
sp_sensitivity_ltz_normal
sp_sensitivity_ltz_uniform01_as_normal
sp_sensitivity_uci_support
sp_uci
spliv
spliv_eval_pattern
spliv_pattern
spliv_sensitivity_path
spliv_tipping_point
```

Registered S3 methods:

```text
plot.plausexog_fit
plot.spliv_fit
plot.spliv_sensitivity_path
print.plausexog_fit
print.spliv_fit
```

## Package tests and check baseline

- `devtools::test()`: 161 passed, 0 failed, 0 warnings, 0 skipped.
- `R CMD build .`: pass; tarball size 57,196 bytes.
- Tarball SHA-256:
  `b99065afd8c9274b85b885b5df53785d94fc3d175b03d205da7a4d0fa6acc0f2`.
- `R CMD check --as-cran spliv_0.2.0.tar.gz`: 0 errors, 0 warnings,
  1 NOTE (maintainer/recent-update incoming-feasibility note; days since last
  update: 2). Examples, tests, vignette rebuild, PDF manual, and HTML manual
  all passed.

## Simulation pilot baseline

Pilot summary checksums:

```text
3dff3c5cecafb478955bf8f77f7df1b12dfc96c8616648f5085959220c4a0b00  simulations/tables/pilot_bpe_summary.csv
05961b16b107b69b8a19c5b8ee7f3d3bc85db884cc067fc9939386a8e8c9c2b5  simulations/tables/pilot_patterned_summary.csv
00de0b21d9892e282286773df68b3b24969d95c834fd8e4c3617653c3e78f9a2  simulations/tables/pilot_subgroup_search_failure.csv
42cbdb223597ad3ba6d1be110f603c943d36772cace8e5f45b07cfa785971c59  simulations/tables/smoke_test_summary.csv
```

Key pilot values before the scaling correction:

- BPE scenarios 1/2/3 eligibility rates: `0.1`, `0.1`, `0.0`.
- BPE scenarios 1/2/3 fit success rates: `0.1`, `0.1`, `0.0`.
- Patterned pilot table: 60 result rows plus header; first scenario coverage
  at delta 0/0.05/0.10/0.15/0.20 was `0.9/1/1/1/1`.
- Subgroup-search false-CI rates for scenarios 1/2/3: `0`, `1`, `0`.
- Smoke summary reported package loaded, baseline and path successful, BPE
  validation ineligible, and no BPE fit.

Full simulation reference checksums before any potential replacement:

```text
a6ac7e077cc50c2dcbf3daf31c9b0d6512d77dfd8802f4f99b287d713752a0e2  results/reference/figures/figure_1_patterned_width_ratios.pdf
5bbd68542c08f08e4a2cce27f223eace7ef3081b47c71321dc65e9b166fc8c13  results/reference/figures/figure_1_patterned_width_ratios.png
3089ec6849873f66d46bc6654c4ca93a523edd82911bf6462a3044ab6819f273  results/reference/figures/figure_2_bpe_coverage.pdf
03e2544bf2a84ac8dc0006ac334b07ac3fc0e7fcac986bf8720a26b16e130670  results/reference/figures/figure_2_bpe_coverage.png
55e6161d556cad6395367b6e6483c691c0c8c0e93450a3391b0f758ade489cdf  results/reference/figures/figure_3_subgroup_search_false_f.pdf
decde1058006f9b70e5624f512e2a097c845829a024db8cf1e547eea50a7eb6c  results/reference/figures/figure_3_subgroup_search_false_f.png
485bc0098cbd263703a044f6341d801b66488bb4a0c95298bc4b60edb8b97965  results/reference/figures/figure_A1_subgroup_search_ci_equivalence.pdf
7d540da173f79bc28b4ab81f8c45bc26db6509d79ef78a56adcc45982ba82778  results/reference/figures/figure_A1_subgroup_search_ci_equivalence.png
2ee9d810d49593740daab7bcf1ed2a910d445101bcd596c5e1fa989819917e16  results/reference/tables/README_paper_sim_outputs.md
6e6346e280ded84db8db8c53c16d1b2e7db61d4becd58e0da89a11eb86211646  results/reference/tables/table_1_patterned_sensitivity_at_true_delta.csv
2dc00d1a88f584f42f994efa2d03ede2613585282fe39f4e37427b45e69b7265  results/reference/tables/table_1_patterned_sensitivity_at_true_delta.tex
89e0139cb7255e930521c4fdeed84b279894649b45482a7b799854e4bb42329f  results/reference/tables/table_2_bpe_design_performance.csv
497b0542948dd5c2b26aec2053abd155d12db12023c2083900d2a745bcae724b  results/reference/tables/table_2_bpe_design_performance.tex
b45cce18d5f6d0ba9ccb40bded1d11e857fb45c259969a4a93fc81ba88040481  results/reference/tables/table_3_subgroup_search_stress.csv
d3bad2c2604667d19a3c69bfc652b3b0a8d97be5d7e71693464dc0509fc85c96  results/reference/tables/table_3_subgroup_search_stress.tex
```

Reference-table dimensions and selected baseline summaries:

- Patterned sensitivity reference: 18 rows; gradient/`pi=0.3`/`theta=0.05`
  correct-pattern coverage `0.969`, width ratio versus uniform `0.6267`.
- BPE reference: 24 rows; zero-gap, `pi=0.3`, inactive share 0.1 had
  eligibility `1`, coverage `0.948`, mean interval width `0.6294`.
- Subgroup-search reference: 72 rows; first row (`pi=0.05`, `K=20`, share
  0.01) false-F, false-CI, and false-equivalence rates were all `1`.

## Koren expected-output baseline

Key results before the standardized BPE-equivalence correction:

- Baseline maize estimate/interval: `184.452590793658`,
  `[74.5001413832304, 294.405040204086]`.
- Baseline wheat estimate/interval: `75.1305028819157`,
  `[29.6703067811097, 120.590698982722]`.
- Uniform UCI tipping-point share: `0.01` for maize and wheat.
- Sparse/bare UCI tipping-point share: `0.03` for maize and wheat.
- Maize BPE: eligible; estimate `182.572581328494`, interval
  `[70.6722936460242, 294.472869010964]`.
- Wheat BPE: ineligible; estimate withheld.
- Old raw-scale equivalence margin: maize `0.000285457998345919`, wheat
  `0.000477597297031356`.
- Old subset raw first-stage CI: maize
  `[-0.0000669108740817623, 0.0000657736389573089]`; wheat
  `[-0.00140970811079699, 0.0000627184561533847]`.

Koren expected-output SHA-256 checksums:

```text
51eef4392981668d29af55c4d704993996e17bd2247f2b79c6a56ad912b347d8  koren_application_summary.csv
2fed144bb05d4e7d519b3525a4b78ed9a04ffbfb9cecdb869123cc01ea290705  koren_baseline_iv.csv
a288df31c88f764273a4ba3fdd4b4a1e7d43c480abfff0ceb9a9f5e2ce9331bf  koren_bpe_diagnostics.csv
879183a409391e6859b03e7cc6dd719daadbd72526a18336cffb1da588ea417e  koren_bpe_diagnostics.pdf
f40bb30e616575d6fad0d3eccde3694cf86d4e696ab602b32a4908d961c70e92  koren_bpe_diagnostics.png
aabb133e72b3a75da64b5113884be75c40e7a8cfaa87df762d1395f6085803a7  koren_ltz_paths.csv
f103cc473d0a69a97ecbbfcd9089a0a393f1568a91cfe38177f6198d87f2f84c  koren_ltz_paths.pdf
713ce75409f79c7f6653ded4b4bddbd232d9ad455ba6e2c56ec31f7ddbabf871  koren_ltz_paths.png
15b1ca303128c9bcb59d8aa622c9d22a9c99347af4dce24f01340a7877e6e45b  koren_ltz_summary.csv
5fd0b210e088b0030bde40acc5a983dda6069b3dcaa3ffe89fb6a6ef7a41e1ea  koren_sensitivity_paths.csv
9c0715e600edfce8a87b6f87279e8a1baa955f62017c11df2dcec52c12a61cc5  koren_sensitivity_paths.pdf
f971c5dbf9a5c02917fc3015e35e96092b8810e1906a5a1f51d7b2df4bb70371  koren_sensitivity_paths.png
5dbe3a6eb68388bc4b5cc8ea3bd262c3dbb42547cfdad2c0b3622aad9d35d0e4  koren_sensitivity_summary.csv
```

No empirical or simulation reference output had been overwritten when this
baseline was written.
