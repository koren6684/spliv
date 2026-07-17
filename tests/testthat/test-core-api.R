test_that("formula parser and instrument names work", {
  d <- make_synth_panel(n_gid = 10, n_t = 10, seed = 11)
  f <- y ~ x + w1 + w2 | z + w1 + w2

  nms <- iv_inst_names(f, d)
  expect_true("z" %in% nms)
  expect_true("(Intercept)" %in% nms)
})

test_that("LTZ and UCI wrappers run on synthetic data", {
  d <- make_synth_panel(n_gid = 20, n_t = 20, seed = 2)
  f <- y ~ x + w1 + w2 | z + w1 + w2

  prior <- conley_prior_ltz(f, d, inst_vary = "z", mean = 0, sd = 0.2)

  ltz <- conley_ltz(
    f,
    d,
    omega = prior$omega,
    mu = prior$mu,
    vcov = "hc1"
  )
  expect_true(is.data.frame(ltz))
  expect_true(all(c("term", "estimate", "conf.low", "conf.high") %in% names(ltz)))

  uci <- conley_uci(
    f,
    d,
    inst = "z",
    gmin = -0.2,
    gmax = 0.2,
    grid = 9,
    vcov = "hc1"
  )
  expect_true(is.data.frame(uci))
  expect_true(all(c("term", "conf.low", "conf.high") %in% names(uci)))
})

test_that("main spliv API supports BPE path", {
  set.seed(3)
  n <- 8000
  z <- rnorm(n)
  u <- rnorm(n)
  e <- rnorm(n)
  inactive_region <- rep(c(TRUE, FALSE), each = n / 2)
  x <- ifelse(inactive_region, 0, 1.2) * z + u
  y <- 1.1 * x + 0.4 * z + e
  d <- data.frame(y = y, x = x, z = z, inactive_region = inactive_region)
  f <- y ~ x - 1 | z - 1
  design <- bpe_design(
    name = "API design",
    subset = ~ inactive_region,
    rationale = "The treatment channel is absent in the inactive region."
  )

  fit <- spliv(
    formula = f,
    data = d,
    method = "bpe",
    vcov = "hc1",
    bpe_design = design,
    bpe_equiv_margin = 0.05,
    bpe_min_n_S = 100
  )

  expect_s3_class(fit, "plausexog_fit")
  expect_true(is.list(fit$prior))
  expect_true(!is.null(fit$bpe_diagnostics$subset_size))
  expect_true("x" %in% fit$estimates$term)
})

test_that("plotting helpers accept fit object", {
  d <- make_synth_panel(n_gid = 15, n_t = 15, seed = 4)
  f <- y ~ x + w1 + w2 | z + w1 + w2

  fit <- spliv(
    formula = f,
    data = d,
    method = "ltz",
    vcov = "hc1",
    prior = conley_prior_ltz(f, d, inst_vary = "z", sd = 0.1)
  )

  sens <- conley_sensitivity_ltz_normal(
    fit,
    term = "x",
    inst_vary = "z",
    delta_grid = seq(0, 0.4, by = 0.2)
  )

  expect_true(is.data.frame(sens))
  expect_true(nrow(sens) == 3)

  rplots_path <- testthat::test_path("Rplots.pdf")
  if (file.exists(rplots_path)) {
    unlink(rplots_path)
  }
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  dev_id <- grDevices::dev.cur()
  on.exit({
    open_devices <- grDevices::dev.list()
    if (!is.null(open_devices) && dev_id %in% open_devices) {
      grDevices::dev.off(dev_id)
    }
    unlink(plot_file)
    if (file.exists(rplots_path)) {
      unlink(rplots_path)
    }
  }, add = TRUE)

  expect_silent(plot_conley_sensitivity(fit))

  path <- spliv_sensitivity_path(
    formula = f,
    data = d,
    method = "uci",
    delta_grid = c(0, 0.1, 0.2),
    vcov = "hc1"
  )
  expect_silent(plot(path, term = "x"))
  expect_false(file.exists(rplots_path))
})

test_that("spliv_pattern accepts formula, function, numeric vector, logical vector, and column-name inputs", {
  d <- make_pattern_sensitivity_data(n = 8, seed = 41)
  numeric_target <- d$treatment_channel_exposure
  logical_target <- as.numeric(d$logical_pattern)

  pattern_formula <- spliv_pattern(
    name = "Formula pattern",
    pattern = ~ treatment_channel_exposure,
    rationale = "Theory implies larger direct effects in higher-exposure areas."
  )
  pattern_function <- spliv_pattern(
    name = "Function pattern",
    pattern = function(data) data$treatment_channel_exposure,
    rationale = "Theory implies larger direct effects in higher-exposure areas.",
    variables_used = "treatment_channel_exposure"
  )
  pattern_numeric <- spliv_pattern(
    name = "Numeric pattern",
    pattern = numeric_target,
    rationale = "Theory implies larger direct effects in higher-exposure areas.",
    normalize = "none"
  )
  pattern_logical <- spliv_pattern(
    name = "Logical pattern",
    pattern = d$logical_pattern,
    rationale = "Theory implies larger direct effects in high-exposure regions.",
    normalize = "none"
  )
  pattern_column <- spliv_pattern(
    name = "Column pattern",
    pattern = "logical_pattern",
    rationale = "Theory implies larger direct effects in high-exposure regions.",
    normalize = "none"
  )

  expect_equal(as.numeric(spliv_eval_pattern(pattern_formula, d)), numeric_target / max(abs(numeric_target)))
  expect_equal(as.numeric(spliv_eval_pattern(pattern_function, d)), numeric_target / max(abs(numeric_target)))
  expect_equal(as.numeric(spliv_eval_pattern(pattern_numeric, d)), numeric_target)
  expect_equal(as.numeric(spliv_eval_pattern(pattern_logical, d)), logical_target)
  expect_equal(as.numeric(spliv_eval_pattern(pattern_column, d)), logical_target)
})

test_that("spliv_eval_pattern returns a numeric vector and ~1 gives a uniform pattern", {
  d <- make_pattern_sensitivity_data(n = 10, seed = 42)
  pat <- spliv_pattern(
    name = "Uniform pattern",
    pattern = ~ 1,
    rationale = "Theory implies a uniform direct effect pattern."
  )

  vals <- spliv_eval_pattern(pat, d)
  expect_type(vals, "double")
  expect_length(vals, nrow(d))
  expect_equal(as.numeric(vals), rep(1, nrow(d)))
})

test_that("logical patterns convert to 0/1 and normalization behaves as requested", {
  d <- make_pattern_sensitivity_data(n = 20, seed = 43)

  pat_logical <- spliv_pattern(
    name = "Logical 0/1",
    pattern = "logical_pattern",
    rationale = "High-exposure regions define the pattern.",
    normalize = "none"
  )
  vals_logical <- spliv_eval_pattern(pat_logical, d)
  expect_true(all(as.numeric(vals_logical) %in% c(0, 1)))

  pat_max_abs <- spliv_pattern(
    name = "Max abs pattern",
    pattern = ~ treatment_channel_exposure,
    rationale = "Theory implies a continuous exposure pattern.",
    normalize = "max_abs"
  )
  vals_max_abs <- spliv_eval_pattern(pat_max_abs, d)
  expect_equal(max(abs(vals_max_abs)), 1, tolerance = 1e-10)

  pat_sd <- spliv_pattern(
    name = "SD pattern",
    pattern = ~ treatment_channel_exposure,
    rationale = "Theory implies a continuous exposure pattern.",
    normalize = "sd"
  )
  vals_sd <- spliv_eval_pattern(pat_sd, d)
  expect_equal(stats::sd(vals_sd), 1, tolerance = 1e-10)

  raw_numeric <- seq_len(nrow(d))
  pat_none <- spliv_pattern(
    name = "Unscaled pattern",
    pattern = raw_numeric,
    rationale = "Theory uses the supplied scale directly.",
    normalize = "none"
  )
  vals_none <- spliv_eval_pattern(pat_none, d)
  expect_equal(as.numeric(vals_none), raw_numeric)
  expect_match(
    attr(vals_none, "spliv_pattern_warnings")[[1]],
    "normalize = 'none'",
    fixed = TRUE
  )
})

test_that("invalid violation patterns fail clearly", {
  d <- make_pattern_sensitivity_data(n = 8, seed = 44)

  expect_error(
    spliv_eval_pattern(
      spliv_pattern(
        name = "Missing variable",
        pattern = ~ missing_exposure,
        rationale = "Invalid pattern."
      ),
      d
    ),
    "Failed to evaluate violation pattern"
  )
  expect_error(
    spliv_eval_pattern(
      spliv_pattern(
        name = "Wrong length",
        pattern = c(1, 2),
        rationale = "Invalid pattern."
      ),
      d
    ),
    "returned length"
  )
  expect_error(
    spliv_eval_pattern(
      spliv_pattern(
        name = "All zero",
        pattern = rep(0, nrow(d)),
        rationale = "Invalid pattern."
      ),
      d
    ),
    "all zero"
  )
  expect_error(
    spliv_eval_pattern(
      spliv_pattern(
        name = "All NA",
        pattern = function(data) rep(NA_real_, nrow(data)),
        rationale = "Invalid pattern.",
        variables_used = "treatment_channel_exposure"
      ),
      d
    ),
    "all `NA`"
  )
  expect_error(
    spliv_eval_pattern(
      spliv_pattern(
        name = "Non-numeric",
        pattern = function(data) rep("a", nrow(data)),
        rationale = "Invalid pattern.",
        variables_used = "treatment_channel_exposure"
      ),
      d
    ),
    "numeric or logical vector"
  )
})

test_that("spliv still works without violation_pattern", {
  d <- make_pattern_sensitivity_data(n = 4000, seed = 45)
  fit <- spliv(
    y ~ x + w | z + w,
    data = d,
    method = "uci",
    vcov = "hc1",
    delta = 0.2
  )

  expect_s3_class(fit, "plausexog_fit")
  expect_true(is.null(fit$violation_pattern))
  expect_equal(fit$grid$parameter, "gamma")
})

test_that("spliv_sensitivity_path works for UCI paths with and without a violation pattern", {
  d <- make_pattern_sensitivity_data(n = 4000, seed = 451)
  f <- y ~ x + w | z + w
  uniform_pattern <- spliv_pattern(
    name = "Uniform direct effect",
    pattern = ~ 1,
    rationale = "Theory implies a uniform direct-effect pattern."
  )
  patterned <- spliv_pattern(
    name = "Exposure pattern",
    pattern = ~ treatment_channel_exposure,
    rationale = paste(
      "The direct effect is expected to be larger where the alternative",
      "channel is more exposed."
    ),
    variables_used = "treatment_channel_exposure",
    pattern_type = "theory_defined",
    normalize = "max_abs"
  )

  default_path <- spliv_sensitivity_path(
    formula = f,
    data = d,
    method = "uci",
    delta_grid = c(0, 0.1, 0.2),
    vcov = "hc1"
  )
  uniform_path <- spliv_sensitivity_path(
    formula = f,
    data = d,
    method = "uci",
    delta_grid = c(0, 0.1, 0.2),
    vcov = "hc1",
    violation_pattern = uniform_pattern
  )
  patterned_path <- spliv_sensitivity_path(
    formula = f,
    data = d,
    method = "uci",
    delta_grid = c(0, 0.1, 0.2),
    vcov = "hc1",
    violation_pattern = patterned
  )

  default_df <- as.data.frame(default_path, stringsAsFactors = FALSE)
  uniform_df <- as.data.frame(uniform_path, stringsAsFactors = FALSE)
  patterned_df <- as.data.frame(patterned_path, stringsAsFactors = FALSE)

  expect_s3_class(default_path, "spliv_sensitivity_path")
  expect_true(inherits(default_path, "data.frame"))
  expect_true(all(
    c(
      "delta", "method", "estimate", "conf_low", "conf_high",
      "contains_zero", "pattern_name", "pattern_type",
      "violation_pattern_used", "scale_instrument", "nobs"
    ) %in% names(default_df)
  ))
  expect_equal(default_df$theta_min, -default_df$delta, tolerance = 1e-12)
  expect_equal(default_df$theta_max, default_df$delta, tolerance = 1e-12)
  expect_equal(attr(default_path, "method"), "uci")
  expect_equal(attr(default_path, "delta_grid"), c(0, 0.1, 0.2))
  expect_true(is.numeric(spliv_tipping_point(default_path)))

  default_x <- default_df[default_df$term == "x", , drop = FALSE]
  uniform_x <- uniform_df[uniform_df$term == "x", , drop = FALSE]
  patterned_x <- patterned_df[patterned_df$term == "x", , drop = FALSE]

  expect_equal(default_x$conf_low, uniform_x$conf_low, tolerance = 1e-10)
  expect_equal(default_x$conf_high, uniform_x$conf_high, tolerance = 1e-10)
  expect_true(all(patterned_df$violation_pattern_used))
  expect_equal(unique(patterned_df$pattern_name), "Exposure pattern")
  expect_equal(unique(patterned_df$pattern_type), "theory_defined")
  expect_gt(max(abs(default_x$conf_low - patterned_x$conf_low)), 1e-6)
  expect_gt(max(abs(default_x$conf_high - patterned_x$conf_high)), 1e-6)
})

test_that("spliv_sensitivity_path works for LTZ paths", {
  d <- make_pattern_sensitivity_data(n = 4000, seed = 452)
  f <- y ~ x + w | z + w
  patterned <- spliv_pattern(
    name = "Exposure pattern",
    pattern = ~ treatment_channel_exposure,
    rationale = "Theory implies larger direct effects in higher-exposure areas.",
    variables_used = "treatment_channel_exposure",
    pattern_type = "theory_defined"
  )

  path <- spliv_sensitivity_path(
    formula = f,
    data = d,
    method = "ltz",
    delta_grid = c(0, 0.05, 0.1),
    vcov = "hc1",
    violation_pattern = patterned
  )
  path_df <- as.data.frame(path, stringsAsFactors = FALSE)

  expect_s3_class(path, "spliv_sensitivity_path")
  expect_true(all(path_df$method == "ltz"))
  expect_true(all(is.na(path_df$theta_min)))
  expect_true(all(is.na(path_df$theta_max)))
  expect_true(all(path_df$violation_pattern_used))
  expect_equal(unique(path_df$pattern_name), "Exposure pattern")
})

test_that("spliv_sensitivity_path computes a tipping point on the supplied grid", {
  d <- make_pattern_sensitivity_data(
    n = 4000,
    beta = 0.05,
    gamma = 0.6,
    seed = 777
  )
  path <- spliv_sensitivity_path(
    formula = y ~ x + w | z + w,
    data = d,
    method = "uci",
    delta_grid = seq(0, 0.5, by = 0.05),
    vcov = "hc1"
  )
  tipping <- spliv_tipping_point(path)

  expect_equal(unname(tipping["x"]), 0.35, tolerance = 1e-12)
  expect_true(is.na(unname(tipping["w"])) || unname(tipping["w"]) >= 0)
  expect_equal(attr(path, "tipping_point")[["x"]], 0.35, tolerance = 1e-12)
})

test_that("spliv_sensitivity_path rejects BPE sensitivity paths", {
  d <- make_pattern_sensitivity_data(n = 3000, seed = 453)

  expect_error(
    spliv_sensitivity_path(
      formula = y ~ x + w | z + w,
      data = d,
      method = "bpe",
      delta_grid = c(0, 0.1)
    ),
    "supports LTZ and UCI only"
  )
})

test_that("uniform violation_pattern matches the default scalar LTZ and UCI paths", {
  d <- make_pattern_sensitivity_data(n = 5000, seed = 46)
  f <- y ~ x + w | z + w
  uniform_pattern <- spliv_pattern(
    name = "Uniform direct effect",
    pattern = ~ 1,
    rationale = "Theory implies a uniform direct effect pattern."
  )

  fit_uci_default <- spliv(f, data = d, method = "uci", vcov = "hc1", delta = 0.2)
  fit_uci_uniform <- spliv(
    f,
    data = d,
    method = "uci",
    vcov = "hc1",
    delta = 0.2,
    violation_pattern = uniform_pattern
  )

  expect_equal(fit_uci_default$estimates$conf.low, fit_uci_uniform$estimates$conf.low, tolerance = 1e-10)
  expect_equal(fit_uci_default$estimates$conf.high, fit_uci_uniform$estimates$conf.high, tolerance = 1e-10)

  fit_ltz_default <- spliv(f, data = d, method = "ltz", vcov = "hc1", delta = 0.2)
  fit_ltz_uniform <- spliv(
    f,
    data = d,
    method = "ltz",
    vcov = "hc1",
    delta = 0.2,
    violation_pattern = uniform_pattern
  )

  expect_equal(fit_ltz_default$estimates$estimate, fit_ltz_uniform$estimates$estimate, tolerance = 1e-10)
  expect_equal(fit_ltz_default$estimates$conf.low, fit_ltz_uniform$estimates$conf.low, tolerance = 1e-10)
  expect_equal(fit_ltz_default$estimates$conf.high, fit_ltz_uniform$estimates$conf.high, tolerance = 1e-10)
})

test_that("non-uniform violation_pattern changes the sensitivity adjustment", {
  d <- make_pattern_sensitivity_data(n = 5000, gamma = 0.6, seed = 47)
  f <- y ~ x + w | z + w
  patterned <- spliv_pattern(
    name = "Exposure pattern",
    pattern = ~ treatment_channel_exposure,
    rationale = paste(
      "The direct effect is expected to be larger where the alternative",
      "channel is more exposed."
    ),
    variables_used = "treatment_channel_exposure",
    pattern_type = "theory_defined",
    normalize = "max_abs"
  )

  fit_default <- spliv(f, data = d, method = "uci", vcov = "hc1", delta = 0.3)
  fit_patterned <- spliv(
    f,
    data = d,
    method = "uci",
    vcov = "hc1",
    delta = 0.3,
    violation_pattern = patterned
  )

  expect_gt(
    max(abs(fit_default$estimates$conf.low - fit_patterned$estimates$conf.low)),
    1e-6
  )
  expect_gt(
    max(abs(fit_default$estimates$conf.high - fit_patterned$estimates$conf.high)),
    1e-6
  )
})

test_that("LTZ and UCI store patterned-sensitivity metadata", {
  d <- make_pattern_sensitivity_data(n = 4000, seed = 48)
  f <- y ~ x + w | z + w
  patterned <- spliv_pattern(
    name = "Theory-defined direct-effect pattern",
    pattern = ~ treatment_channel_exposure,
    rationale = paste(
      "The direct effect is expected to be larger where the alternative",
      "channel is more exposed."
    ),
    variables_used = "treatment_channel_exposure",
    pattern_type = "theory_defined",
    normalize = "max_abs"
  )

  fit_ltz <- spliv(
    f,
    data = d,
    method = "ltz",
    vcov = "hc1",
    delta = 0.1,
    violation_pattern = patterned
  )
  fit_uci <- spliv(
    f,
    data = d,
    method = "uci",
    vcov = "hc1",
    delta = 0.1,
    violation_pattern = patterned
  )

  expect_equal(fit_ltz$violation_pattern$name, "Theory-defined direct-effect pattern")
  expect_equal(fit_ltz$violation_pattern$normalize, "max_abs")
  expect_true(is.list(fit_ltz$violation_pattern$raw_pattern_summary))
  expect_true(is.list(fit_ltz$violation_pattern$normalized_pattern_summary))
  expect_true(is.list(fit_ltz$diag$violation_pattern))
  expect_true(fit_ltz$violation_pattern$residual_sd_scaling_used)

  expect_equal(fit_uci$grid$gmin, -0.1)
  expect_equal(fit_uci$grid$gmax, 0.1)
  expect_equal(fit_uci$grid$parameter, "theta")
  expect_equal(fit_uci$violation_pattern$instrument, "z")
  expect_equal(fit_uci$violation_pattern$endogenous_treatment, "x")
})

test_that("patterned sensitivity errors clearly for unsupported multiple targets", {
  d_multi_z <- make_pattern_sensitivity_data(n = 3000, seed = 49)
  d_multi_z$z2 <- d_multi_z$z + rnorm(nrow(d_multi_z))
  patterned <- spliv_pattern(
    name = "Exposure pattern",
    pattern = ~ treatment_channel_exposure,
    rationale = "Theory implies larger direct effects in higher-exposure areas."
  )

  expect_error(
    spliv(
      y ~ x + w | z + z2 + w,
      data = d_multi_z,
      method = "uci",
      vcov = "hc1",
      delta = 0.2,
      violation_pattern = patterned
    ),
    "one endogenous treatment and one instrument"
  )

  d_multi_x <- make_pattern_sensitivity_data(n = 3000, seed = 50)
  d_multi_x$x2 <- 0.7 * d_multi_x$z + 0.3 * d_multi_x$w + rnorm(nrow(d_multi_x))
  expect_error(
    spliv(
      y ~ x + x2 + w | z + w,
      data = d_multi_x,
      method = "ltz",
      vcov = "hc1",
      delta = 0.2,
      violation_pattern = patterned
    ),
    "one endogenous treatment and one instrument"
  )
})

test_that("confirmatory BPE intentionally rejects violation_pattern", {
  d <- make_bpe_redesign_data(n = 4000, seed = 51)
  design <- bpe_design(
    name = "Inactive subset",
    subset = ~ inactive_region,
    rationale = "The treatment channel is absent in the inactive region."
  )
  patterned <- spliv_pattern(
    name = "Uniform pattern",
    pattern = ~ 1,
    rationale = "A uniform direct-effect pattern."
  )

  expect_error(
    spliv(
      y ~ x - 1 | z - 1,
      data = d,
      method = "bpe",
      bpe_design = design,
      bpe_equiv_margin = 0.05,
      violation_pattern = patterned
    ),
    "not currently supported with confirmatory BPE"
  )
})

test_that("bpe_design accepts formula, function, logical vector, and column-name inputs", {
  d <- make_bpe_redesign_data(n = 8, seed = 11)
  target <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)

  design_formula <- bpe_design(
    name = "Formula design",
    subset = ~ theoretical_condition == 1,
    rationale = "Theory identifies an inactive subset."
  )
  design_function <- bpe_design(
    name = "Function design",
    subset = function(data) data$inactive_region,
    rationale = "Theory identifies an inactive subset.",
    variables_used = "inactive_region"
  )
  design_logical <- bpe_design(
    name = "Logical design",
    subset = target,
    rationale = "Theory identifies an inactive subset."
  )
  design_column <- bpe_design(
    name = "Column design",
    subset = "logical_subset",
    rationale = "Theory identifies an inactive subset."
  )

  expect_equal(as.vector(bpe_eval_subset(design_formula, d)), as.vector(target))
  expect_equal(as.vector(bpe_eval_subset(design_function, d)), as.vector(target))
  expect_equal(as.vector(bpe_eval_subset(design_logical, d)), as.vector(target))
  expect_equal(as.vector(bpe_eval_subset(design_column, d)), as.vector(target))
})

test_that("invalid BPE subsets fail clearly", {
  d <- make_bpe_redesign_data(n = 8, seed = 12)

  expect_error(
    bpe_eval_subset(
      bpe_design("wrong length", c(TRUE, FALSE), "rationale"),
      d
    ),
    "returned length"
  )
  expect_error(
    bpe_eval_subset(
      bpe_design("all true", rep(TRUE, nrow(d)), "rationale"),
      d
    ),
    "selected all observations"
  )
  expect_error(
    bpe_eval_subset(
      bpe_design("all false", rep(FALSE, nrow(d)), "rationale"),
      d
    ),
    "selected no observations"
  )
  expect_error(
    bpe_eval_subset(
      bpe_design("missing var", ~ missing_region, "rationale"),
      d
    ),
    "Failed to evaluate"
  )
  expect_error(
    bpe_eval_subset(
      bpe_design(
        "non-logical",
        function(data) seq_len(nrow(data)),
        "rationale",
        variables_used = "inactive_region"
      ),
      d
    ),
    "must evaluate to a logical vector"
  )
})

test_that("confirmatory BPE fails when no design object is supplied", {
  d <- make_bpe_redesign_data(n = 4000, seed = 13)

  expect_error(
    spliv(
      y ~ x - 1 | z - 1,
      data = d,
      method = "bpe",
      bpe_equiv_margin = 0.05
    ),
    "requires a researcher-supplied `bpe_design\\(\\)` object"
  )
})

test_that("confirmatory BPE errors when multiple subset sources are supplied", {
  d <- make_bpe_redesign_data(n = 4000, seed = 131)
  design_1 <- bpe_design(
    name = "Source one",
    subset = ~ inactive_region,
    rationale = "Theory identifies the inactive subset."
  )
  design_2 <- bpe_design(
    name = "Source two",
    subset = ~ theoretical_condition == 1,
    rationale = "A second declared source should trigger an error."
  )

  expect_error(
    spliv(
      y ~ x - 1 | z - 1,
      data = d,
      method = "bpe",
      bpe_design = design_1,
      bpe_spec = list(design = design_2),
      bpe_equiv_margin = 0.05
    ),
    "exactly one declared design/subset source"
  )
})

test_that("raw bpe_spec subset requires explicit confirmatory metadata", {
  d <- make_bpe_redesign_data(n = 5000, seed = 132)

  expect_error(
    spliv(
      y ~ x - 1 | z - 1,
      data = d,
      method = "bpe",
      bpe_spec = list(
        subset = ~ treatment_channel_absent,
        rationale = "Theory identifies the inactive subset."
      ),
      bpe_equiv_margin = 0.05
    ),
    "explicit `pre_specified = TRUE`"
  )

  expect_error(
    spliv(
      y ~ x - 1 | z - 1,
      data = d,
      method = "bpe",
      bpe_spec = list(
        subset = ~ treatment_channel_absent,
        rationale = "Theory identifies the inactive subset.",
        pre_specified = FALSE
      ),
      bpe_equiv_margin = 0.05
    ),
    "requires `pre_specified = TRUE`"
  )

  fit <- spliv(
    y ~ x - 1 | z - 1,
    data = d,
    method = "bpe",
    bpe_spec = list(
      subset = ~ treatment_channel_absent,
      rationale = "Theory identifies the inactive subset.",
      pre_specified = TRUE,
      variables_used = "treatment_channel_absent"
    ),
    bpe_equiv_margin = 0.05,
    bpe_min_n_S = 100
  )

  expect_s3_class(fit, "plausexog_fit")
})

test_that("confirmatory BPE fails when pre_specified is FALSE", {
  d <- make_bpe_redesign_data(n = 4000, seed = 14)
  design <- bpe_design(
    name = "Exploratory subset",
    subset = ~ inactive_region,
    rationale = "This subset was found after looking at the data.",
    pre_specified = FALSE
  )

  expect_error(
    spliv(
      y ~ x - 1 | z - 1,
      data = d,
      method = "bpe",
      bpe_design = design,
      bpe_equiv_margin = 0.05
    ),
    "pre_specified = TRUE"
  )
})

test_that("confirmatory BPE fails when rationale is missing", {
  d <- make_bpe_redesign_data(n = 4000, seed = 15)
  design <- bpe_design(
    name = "Missing rationale",
    subset = ~ inactive_region,
    rationale = ""
  )

  expect_error(
    spliv(
      y ~ x - 1 | z - 1,
      data = d,
      method = "bpe",
      bpe_design = design,
      bpe_equiv_margin = 0.05
    ),
    "non-empty `rationale`"
  )
})

test_that("confirmatory BPE errors when multiple endogenous regressors are present", {
  set.seed(151)
  n <- 5000
  z <- rnorm(n)
  inactive_region <- rep(c(TRUE, FALSE), each = n / 2)
  x1 <- ifelse(inactive_region, 0.02, 0.9) * z + rnorm(n)
  x2 <- ifelse(inactive_region, 0.01, 0.7) * z + rnorm(n)
  y <- 1.1 * x1 - 0.3 * x2 + 0.4 * z + rnorm(n)
  d <- data.frame(y = y, x1 = x1, x2 = x2, z = z, inactive_region = inactive_region)
  design <- bpe_design(
    name = "Multiple endogenous regressors",
    subset = ~ inactive_region,
    rationale = "Theory identifies the inactive subset."
  )

  expect_error(
    bpe_validate_design(
      y ~ x1 + x2 - 1 | z - 1,
      data = d,
      design = design,
      bpe_equiv_margin = 0.05
    ),
    "supports one endogenous treatment"
  )
})

test_that("first-stage coefficient extraction works with and without intercepts and fixed effects", {
  d_no_intercept <- make_bpe_redesign_data(n = 8000, pi_S = 0.08, pi_notS = 1.0, seed = 16)
  design_no_intercept <- bpe_design(
    name = "No intercept",
    subset = ~ inactive_region,
    rationale = "The instrument is inactive in this design-defined subset."
  )

  val_no_intercept <- bpe_validate_design(
    y ~ x - 1 | z - 1,
    data = d_no_intercept,
    design = design_no_intercept,
    bpe_equiv_margin = 0.2,
    bpe_min_n_S = 100
  )

  expect_true("z" %in% names(val_no_intercept$first_stage_coefficient))
  expect_equal(unname(val_no_intercept$first_stage_coefficient["z"]), 0.08, tolerance = 0.03)

  skip_if_not_installed("fixest")
  set.seed(17)
  n_id <- 80
  TT <- 8
  n <- n_id * TT
  id <- factor(rep(seq_len(n_id), each = TT))
  tt <- factor(rep(seq_len(TT), times = n_id))
  inactive_region <- rep(rep(c(TRUE, FALSE), each = TT / 2), length.out = n)
  z <- rnorm(n)
  fe_i <- rnorm(n_id)[id]
  fe_t <- rnorm(TT)[tt]
  x <- ifelse(inactive_region, 0.12, 1.1) * z + fe_i + fe_t + rnorm(n)
  y <- 1.1 * x + 0.3 * z + fe_i + fe_t + rnorm(n)
  d_fe <- data.frame(y = y, x = x, z = z, id = id, tt = tt, inactive_region = inactive_region)
  design_fe <- bpe_design(
    name = "FE subset",
    subset = ~ inactive_region,
    rationale = "The channel is absent in the inactive region cells."
  )

  val_fe <- bpe_validate_design(
    y ~ x | z,
    data = d_fe,
    design = design_fe,
    fe = ~ id + tt,
    bpe_equiv_margin = 0.3,
    bpe_min_n_S = 100
  )

  expect_true("z" %in% names(val_fe$first_stage_coefficient))
  expect_lt(abs(unname(val_fe$first_stage_coefficient["z"]) - 0.12), 0.05)
})

test_that("first-stage F-statistic is reported but does not determine eligibility", {
  set.seed(18)
  n <- 20000
  z <- rnorm(n)
  inactive_region <- rep(c(TRUE, FALSE), each = n / 2)
  x <- ifelse(inactive_region, 0.01, 1.1) * z + rnorm(n, sd = 0.1)
  y <- 1.2 * x + 0.4 * z + rnorm(n, sd = 0.1)
  d <- data.frame(y = y, x = x, z = z, inactive_region = inactive_region)
  design <- bpe_design(
    name = "High-power inactive subset",
    subset = ~ inactive_region,
    rationale = "The treatment channel is absent in the inactive region."
  )

  val <- bpe_validate_design(
    y ~ x - 1 | z - 1,
    data = d,
    design = design,
    bpe_equiv_margin = 0.02,
    bpe_min_n_S = 1000
  )

  expect_gt(unname(val$first_stage_f_statistic["z"]), 5)
  expect_true(val$equivalence_passed)
  expect_true(val$eligibility_passed)
})

test_that("equivalence check determines eligibility", {
  d <- make_bpe_redesign_data(
    n = 10000,
    pi_S = 0.03,
    pi_notS = 1.1,
    seed = 19
  )
  design <- bpe_design(
    name = "Tight equivalence subset",
    subset = ~ inactive_region,
    rationale = "The treatment channel is absent in the inactive region."
  )

  val <- bpe_validate_design(
    y ~ x - 1 | z - 1,
    data = d,
    design = design,
    bpe_equiv_margin = 0.005,
    bpe_min_n_S = 1000
  )

  expect_false(val$equivalence_passed)
  expect_false(val$eligibility_passed)
})

test_that("cluster-count check works when clusters are supplied", {
  set.seed(20)
  gid <- factor(rep(seq_len(12), each = 100))
  inactive_region <- gid %in% levels(gid)[1:5]
  z <- rnorm(length(gid))
  x <- ifelse(inactive_region, 0, 1.1) * z + rnorm(length(gid))
  y <- 1.2 * x + 0.4 * z + rnorm(length(gid))
  d <- data.frame(y = y, x = x, z = z, gid = gid, inactive_region = inactive_region)
  design <- bpe_design(
    name = "Cluster-limited subset",
    subset = ~ inactive_region,
    rationale = "The treatment channel is absent in the inactive subset."
  )

  val <- bpe_validate_design(
    y ~ x - 1 | z - 1,
    data = d,
    design = design,
    vcov = "cluster",
    cluster = ~ gid,
    bpe_equiv_margin = 0.05,
    bpe_min_n_S = 100,
    bpe_min_clusters_S = 6
  )

  expect_equal(val$G_S, 5)
  expect_false(val$eligibility_checks$minimum_clusters)
  expect_false(val$eligibility_passed)
})

test_that("reduced-form covariance matrix is stored and propagated", {
  d <- make_bpe_redesign_data(n = 10000, pi_S = 0, pi_notS = 1.0, seed = 21)
  design <- bpe_design(
    name = "Stored covariance subset",
    subset = ~ inactive_region,
    rationale = "The treatment channel is absent in the inactive region.",
    transportability_rationale = "The direct effect is assumed informative for the target sample."
  )

  fit <- spliv(
    y ~ x - 1 | z - 1,
    data = d,
    method = "bpe",
    bpe_design = design,
    bpe_kappa = 1,
    bpe_transport = "sampling",
    bpe_equiv_margin = 0.05,
    bpe_min_n_S = 1000
  )

  z_idx <- match("z", colnames(fit$internals$Z))
  expect_true(is.matrix(fit$bpe_diagnostics$reduced_form_direct_effect_cov))
  expect_equal(dim(fit$bpe_diagnostics$reduced_form_direct_effect_cov), c(1, 1))
  expect_equal(fit$bpe_diagnostics$prior_Omega_sub, fit$bpe_diagnostics$reduced_form_direct_effect_cov)
  expect_equal(
    fit$Omega_used[z_idx, z_idx, drop = FALSE],
    fit$bpe_diagnostics$prior_Omega_sub,
    tolerance = 1e-10
  )
})

test_that("confirmatory BPE works with fixed effects and clustered covariance", {
  set.seed(4201)
  n_unit <- 80
  n_period <- 5
  n <- n_unit * n_period
  unit <- factor(rep(seq_len(n_unit), each = n_period))
  period <- factor(rep(seq_len(n_period), times = n_unit))
  inactive_region <- rep(seq_len(n_unit) <= 35, each = n_period)
  z <- stats::rnorm(n)
  unit_x <- stats::rnorm(n_unit)[as.integer(unit)]
  period_x <- stats::rnorm(n_period)[as.integer(period)]
  unit_y <- stats::rnorm(n_unit)[as.integer(unit)]
  period_y <- stats::rnorm(n_period)[as.integer(period)]
  x <- ifelse(inactive_region, 0, 1) * z + unit_x + period_x + stats::rnorm(n, sd = 0.35)
  y <- 1.4 * x + 0.15 * z + unit_y + period_y + stats::rnorm(n, sd = 0.35)
  d <- data.frame(
    y = y,
    x = x,
    z = z,
    unit = unit,
    period = period,
    inactive_region = inactive_region
  )
  design <- bpe_design(
    name = "Koren-like inactive subset",
    subset = ~ inactive_region,
    rationale = "The treatment channel is absent in this pre-specified inactive subset.",
    variables_used = "inactive_region",
    pre_specified = TRUE
  )

  fit <- spliv(
    y ~ x | z,
    data = d,
    fe = ~ unit + period,
    vcov = "cluster",
    cluster = ~ unit,
    method = "bpe",
    bpe_design = design,
    bpe_equiv_margin = 0.25,
    bpe_min_n_S = 100,
    bpe_min_clusters_S = 20,
    scale_instrument = "residual_sd"
  )

  expect_s3_class(fit, "plausexog_fit")
  expect_true(isTRUE(fit$bpe_diagnostics$eligibility_passed))
  expect_equal(length(fit$mu_used), ncol(fit$internals$Z))
  expect_equal(dim(fit$Omega_used), c(ncol(fit$internals$Z), ncol(fit$internals$Z)))
  expect_true(all(is.finite(fit$estimates$estimate)))
  expect_true(all(is.finite(fit$estimates$conf.low)))
  expect_true(all(is.finite(fit$estimates$conf.high)))
})

test_that("standardized first-stage diagnostics are stored when available", {
  d <- make_bpe_redesign_data(n = 10000, pi_S = 0.08, pi_notS = 1.0, seed = 211)
  design <- bpe_design(
    name = "Scaled first stage subset",
    subset = ~ inactive_region,
    rationale = "The treatment channel is absent in the inactive region."
  )

  val <- bpe_validate_design(
    y ~ x - 1 | z - 1,
    data = d,
    design = design,
    bpe_equiv_margin = 0.2,
    bpe_min_n_S = 1000
  )

  expect_true(is.finite(unname(val$residualized_instrument_sd_S["z"])))
  expect_true(is.finite(unname(val$residualized_treatment_sd_S["x"])))
  expect_true(is.finite(unname(val$first_stage_effect_one_residual_sd_Z["z"])))
  expect_true(is.finite(unname(val$standardized_first_stage_effect["z"])))
  expect_equal(
    unname(val$first_stage_effect_one_residual_sd_Z["z"]),
    unname(val$first_stage_coefficient["z"] * val$residualized_instrument_sd_S["z"]),
    tolerance = 1e-10
  )
  expect_equal(
    unname(val$standardized_first_stage_effect["z"]),
    unname(val$first_stage_effect_one_residual_sd_Z["z"] / val$residualized_treatment_sd_S["x"]),
    tolerance = 1e-10
  )
})

test_that("scale_instrument = 'residual_sd' stores the residualized instrument SD", {
  d <- make_bpe_redesign_data(n = 10000, pi_S = 0, pi_notS = 1.0, seed = 22)
  design <- bpe_design(
    name = "Scaled subset",
    subset = ~ inactive_region,
    rationale = "The treatment channel is absent in the inactive region."
  )

  fit <- spliv(
    y ~ x - 1 | z - 1,
    data = d,
    method = "bpe",
    bpe_design = design,
    bpe_equiv_margin = 0.05,
    scale_instrument = "residual_sd",
    bpe_min_n_S = 1000
  )

  expect_true(is.finite(fit$residualized_instrument_sd["z"]))
  expect_true(is.finite(fit$bpe_diagnostics$residualized_instrument_sd["z"]))
  expect_equal(
    unname(fit$residualized_instrument_sd["z"]),
    unname(fit$bpe_diagnostics$residualized_instrument_sd["z"]),
    tolerance = 1e-10
  )
})

test_that("UCI delta defaults correspond to [-delta, +delta]", {
  d <- make_synth_panel(n_gid = 20, n_t = 20, seed = 23)
  f <- y ~ x + w1 + w2 | z + w1 + w2

  default_path <- conley_sensitivity_uci_support(
    formula = f,
    data = d,
    term = "x",
    inst_vary = "z",
    delta_grid = c(0.2, 0.4),
    scale_instrument = "none"
  )
  explicit_path <- conley_sensitivity_uci_support(
    formula = f,
    data = d,
    term = "x",
    inst_vary = "z",
    delta_grid = c(0.2, 0.4),
    gmin_fun = function(delta) -delta,
    gmax_fun = function(delta) delta,
    scale_instrument = "none"
  )

  expect_equal(default_path$conf.low, explicit_path$conf.low, tolerance = 1e-12)
  expect_equal(default_path$conf.high, explicit_path$conf.high, tolerance = 1e-12)
})
