test_that("plausexog wrapper is deprecated and matches spliv outputs", {
  d <- make_synth_panel(n_gid = 40, n_t = 20, seed = 20260223)
  f <- y ~ x + w1 + w2 | z + w1 + w2

  expect_warning(
    plausexog(
      formula = f,
      data = d,
      method = "uci",
      vcov = "hc1",
      grid = list(inst = "z", gmin = -0.2, gmax = 0.2, steps = 5)
    ),
    "deprecated"
  )

  prior <- conley_prior_ltz(f, d, inst_vary = "z", mean = 0, sd = 0.2)
  fit_ltz_new <- spliv(
    formula = f,
    data = d,
    method = "ltz",
    vcov = "hc1",
    prior = prior
  )
  fit_ltz_old <- suppressWarnings(
    plausexog(
      formula = f,
      data = d,
      method = "ltz",
      vcov = "hc1",
      prior = prior
    )
  )
  expect_equal(fit_ltz_new$estimates$estimate, fit_ltz_old$estimates$estimate, tolerance = 1e-12)
  expect_equal(fit_ltz_new$estimates$conf.low, fit_ltz_old$estimates$conf.low, tolerance = 1e-12)
  expect_equal(fit_ltz_new$estimates$conf.high, fit_ltz_old$estimates$conf.high, tolerance = 1e-12)

  fit_uci_new <- spliv(
    formula = f,
    data = d,
    method = "uci",
    vcov = "hc1",
    grid = list(inst = "z", gmin = -0.2, gmax = 0.2, steps = 5)
  )
  fit_uci_old <- suppressWarnings(
    plausexog(
      formula = f,
      data = d,
      method = "uci",
      vcov = "hc1",
      grid = list(inst = "z", gmin = -0.2, gmax = 0.2, steps = 5)
    )
  )
  expect_equal(fit_uci_new$estimates$conf.low, fit_uci_old$estimates$conf.low, tolerance = 1e-12)
  expect_equal(fit_uci_new$estimates$conf.high, fit_uci_old$estimates$conf.high, tolerance = 1e-12)

  set.seed(20260224)
  n <- 6000
  z <- rnorm(n)
  u <- rnorm(n)
  e <- rnorm(n)
  inactive_region <- rep(c(TRUE, FALSE), each = n / 2)
  d_bpe <- data.frame(
    y = 1.1 * (ifelse(inactive_region, 0, 1.0) * z + u) + 0.4 * z + e,
    x = ifelse(inactive_region, 0, 1.0) * z + u,
    z = z,
    inactive_region = inactive_region
  )
  f_bpe <- y ~ x - 1 | z - 1
  design <- bpe_design(
    name = "Deprecated-wrapper design",
    subset = ~ inactive_region,
    rationale = "The treatment channel is absent in the inactive region."
  )

  fit_bpe_new <- spliv(
    formula = f_bpe,
    data = d_bpe,
    method = "bpe",
    vcov = "hc1",
    bpe_design = design,
    bpe_equiv_margin = 0.05,
    bpe_min_n_S = 100
  )
  fit_bpe_old <- suppressWarnings(
    plausexog(
      formula = f_bpe,
      data = d_bpe,
      method = "bpe",
      vcov = "hc1",
      bpe_design = design,
      bpe_equiv_margin = 0.05,
      bpe_min_n_S = 100
    )
  )
  expect_equal(fit_bpe_new$estimates$estimate, fit_bpe_old$estimates$estimate, tolerance = 1e-12)
  expect_equal(fit_bpe_new$estimates$conf.low, fit_bpe_old$estimates$conf.low, tolerance = 1e-12)
  expect_equal(fit_bpe_new$estimates$conf.high, fit_bpe_old$estimates$conf.high, tolerance = 1e-12)
})

test_that("bpe_explore_subsets has no application-specific default rules", {
  expect_null(formals(bpe_explore_subsets)$rules)

  body_text <- paste(deparse(body(bpe_explore_subsets)), collapse = " ")
  expect_false(grepl("sparsebare|crop_area|irrigation|maize|wheat|rainfall", body_text, ignore.case = TRUE))

  d <- make_bpe_redesign_data(n = 2000, seed = 301)
  expect_error(
    suppressWarnings(
      bpe_explore_subsets(
        data = d,
        spec = y ~ x - 1 | z - 1
      )
    ),
    "`rules` must supply"
  )
})

test_that("bpe_explore_subsets warns that it is exploratory only", {
  d <- make_bpe_redesign_data(n = 2000, seed = 302)

  expect_warning(
    out <- bpe_explore_subsets(
      data = d,
      spec = y ~ x - 1 | z - 1,
      rules = list(
        list(name = "formula_rule", subset = ~ inactive_region),
        list(name = "column_rule", subset = "logical_subset"),
        list(
          name = "function_rule",
          subset = function(data) data$treatment_channel_absent
        )
      )
    ),
    "This function is exploratory"
  )

  expect_s3_class(out, "data.frame")
  expect_equal(out$rule, c("formula_rule", "column_rule", "function_rule"))
})

test_that("estimate_gamma_zero_first_stage is deprecated as legacy exploratory", {
  d <- make_bpe_redesign_data(n = 4000, pi_S = 0, pi_notS = 1, seed = 303)

  expect_warning(
    gamma_hat <- estimate_gamma_zero_first_stage(
      data = d,
      y_name = "y",
      z_names = "z",
      subset = ~ inactive_region
    ),
    "legacy/exploratory helper"
  )

  expect_true(all(c("mu_hat", "omega_hat", "diagnostics") %in% names(gamma_hat)))
  expect_true(is.list(gamma_hat$diagnostics))
  expect_true(is.finite(gamma_hat$diagnostics$n_subset))
})

test_that("pattern defaults and README remain generic", {
  pattern_body <- paste(deparse(body(spliv_pattern)), collapse = " ")
  expect_false(grepl("crop|sparse_bare|crop_area|irrigation|maize|wheat|rainfall|koren|ritter", pattern_body, ignore.case = TRUE))

  pkg_desc <- utils::packageDescription("spliv")
  desc_text <- paste(pkg_desc$Title, pkg_desc$Description, collapse = "\n")
  expect_false(grepl("crop|sparse_bare|crop_area|irrigation|maize|wheat|rainfall|drought|koren|ritter", desc_text, ignore.case = TRUE))
  expect_false(grepl("spatiotemporal", desc_text, ignore.case = TRUE))
})
