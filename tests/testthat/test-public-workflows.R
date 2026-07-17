test_that("README-style synthetic path workflow runs", {
  d <- make_pattern_sensitivity_data(n = 1200, seed = 9001)
  f <- y ~ x + w | z + w
  fit <- spliv(f, d, method = "uci", delta = 0.05, vcov = "hc1")
  expect_s3_class(fit, "plausexog_fit")
  pat <- spliv_pattern(
    "Exposure", ~ treatment_channel_exposure,
    rationale = "Synthetic exposure pattern.", variables_used = "treatment_channel_exposure"
  )
  path <- spliv_sensitivity_path(
    f, d, method = "uci", delta_grid = c(0, 0.05),
    violation_pattern = pat, vcov = "hc1"
  )
  expect_s3_class(path, "spliv_sensitivity_path")
  expect_true(all(c("estimate", "conf_low", "conf_high") %in% names(path)))
  expect_type(spliv_tipping_point(path), "double")
})

test_that("confirmatory BPE design workflow is callable on synthetic data", {
  set.seed(9002)
  n <- 1200
  z <- rnorm(n)
  inactive <- rep(c(TRUE, FALSE), each = n / 2)
  x <- ifelse(inactive, 0, 1) * z + rnorm(n)
  y <- 1.1 * x + 0.2 * z + rnorm(n)
  d <- data.frame(y, x, z, inactive)
  design <- bpe_design(
    "Synthetic inactive subset", ~ inactive,
    rationale = "The treatment channel is absent in the inactive subset.",
    variables_used = "inactive", pre_specified = TRUE
  )
  val <- bpe_validate_design(
    y ~ x | z, d, design, vcov = "hc1",
    bpe_min_n_S = 100, bpe_equiv_margin = 1
  )
  expect_s3_class(val, "spliv_bpe_validation")
  expect_true(is.list(val))
})
