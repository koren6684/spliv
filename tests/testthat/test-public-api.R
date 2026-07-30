test_that("the exported API is exactly the approved interface", {
  canonical <- c(
    "spliv",
    "spliv_pattern",
    "spliv_eval_pattern",
    "spliv_sensitivity_path",
    "spliv_tipping_point",
    "bpe_design",
    "bpe_eval_subset",
    "bpe_validate_design",
    "bpe_explore_subsets"
  )
  advanced <- c(
    "sp_ltz",
    "sp_uci",
    "sp_prior_ltz",
    "sp_sensitivity_ltz_normal",
    "sp_sensitivity_uci_support",
    "sp_sensitivity_ltz_uniform01_as_normal",
    "plot_sp_sensitivity"
  )
  expected_exports <- sort(c(canonical, advanced))
  exports <- sort(getNamespaceExports("spliv"))

  expect_true(all(canonical %in% exports))
  expect_identical(exports, expected_exports)
})
