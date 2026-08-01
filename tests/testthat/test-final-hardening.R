test_that("baseline IV agrees with fixest across FE and covariance choices", {
  skip_if_not_installed("fixest")
  d <- make_synth_panel(n_gid = 30, n_t = 20, seed = 3101)
  f <- y ~ x + w1 + w2 | z + w1 + w2
  ssc <- fixest::ssc(
    K.adj = TRUE, K.fixef = "none", G.adj = TRUE,
    G.df = "min", t.df = "min"
  )

  for (with_fe in c(FALSE, TRUE)) {
    for (vcov_i in c("iid", "hc1", "cluster")) {
      fe_i <- if (with_fe) ~ gid + ym else NULL
      cluster_i <- if (identical(vcov_i, "cluster")) ~ gid else NULL
      fit <- spliv(
        f, d, fe = fe_i, method = "uci", delta = 0,
        vcov = vcov_i, cluster = cluster_i
      )
      ref_vcov <- switch(vcov_i, iid = "iid", hc1 = "hetero", cluster = ~ gid)
      ref <- if (with_fe) {
        fixest::feols(
          y ~ w1 + w2 | gid + ym | x ~ z,
          data = d, vcov = ref_vcov, ssc = ssc
        )
      } else {
        fixest::feols(
          y ~ w1 + w2 | 0 | x ~ z,
          data = d, vcov = ref_vcov, ssc = ssc
        )
      }

      ref_names <- if (with_fe) c("fit_x", "w1", "w2") else c("(Intercept)", "fit_x", "w1", "w2")
      ref_idx <- match(ref_names, names(stats::coef(ref)))
      fit_estimate <- (fit$estimates$conf.low + fit$estimates$conf.high) / 2
      fit_se <- (fit$estimates$conf.high - fit$estimates$conf.low) /
        (2 * stats::qnorm(0.975))

      expect_equal(fit_estimate, unname(stats::coef(ref)[ref_idx]), tolerance = 1e-9)
      expect_equal(fit_se, unname(sqrt(diag(stats::vcov(ref)))[ref_idx]), tolerance = 1e-9)
    }
  }
})

test_that("manual adjusted-outcome refits reproduce point and bounded UCI intervals", {
  skip_if_not_installed("fixest")
  set.seed(3102)
  n <- 700
  z <- rnorm(n)
  w <- rnorm(n)
  x <- 0.8 * z + 0.3 * w + rnorm(n)
  y <- 1.1 * x + 0.25 * w + 0.12 * z + rnorm(n)
  d <- data.frame(y, x, z, w)
  theta <- 0.12
  theta_grid <- seq(-0.2, 0.2, length.out = 5)
  ssc <- fixest::ssc(K.adj = TRUE, K.fixef = "none")

  manual_interval <- function(theta_i) {
    d$adjusted_y <- d$y - theta_i * d$z
    ref <- fixest::feols(
      adjusted_y ~ w | 0 | x ~ z,
      data = d, vcov = "hetero", ssc = ssc
    )
    estimate <- unname(stats::coef(ref)[["fit_x"]])
    se <- unname(sqrt(diag(stats::vcov(ref)))[["fit_x"]])
    c(lower = estimate - stats::qnorm(0.975) * se,
      upper = estimate + stats::qnorm(0.975) * se)
  }

  point <- spliv(
    y ~ x + w | z + w, d, method = "uci", vcov = "hc1",
    scale_instrument = "none",
    grid = list(gmin = theta, gmax = theta, steps = 1)
  )
  point_ref <- manual_interval(theta)
  point_row <- point$estimates[point$estimates$term == "x", ]
  expect_equal(point_row$conf.low, unname(point_ref[["lower"]]), tolerance = 1e-9)
  expect_equal(point_row$conf.high, unname(point_ref[["upper"]]), tolerance = 1e-9)

  bounded <- spliv(
    y ~ x + w | z + w, d, method = "uci", vcov = "hc1",
    scale_instrument = "none",
    grid = list(gmin = min(theta_grid), gmax = max(theta_grid), steps = length(theta_grid))
  )
  manual_grid <- vapply(theta_grid, manual_interval, numeric(2))
  bounded_row <- bounded$estimates[bounded$estimates$term == "x", ]
  expect_equal(bounded_row$conf.low, min(manual_grid["lower", ]), tolerance = 1e-9)
  expect_equal(bounded_row$conf.high, max(manual_grid["upper", ]), tolerance = 1e-9)
})

test_that("standardized BPE diagnostics and final estimates are invariant to instrument units", {
  d <- make_bpe_redesign_data(n = 6000, pi_S = 0, pi_notS = 1, seed = 3103)
  design <- bpe_design(
    "Inactive subset", ~ inactive_region,
    rationale = "The treatment channel is absent in this pre-specified subset.",
    transportability_rationale = "The direct-effect mechanism is assumed to apply to the target sample."
  )

  run_scaled <- function(multiplier) {
    scaled <- d
    scaled$z_scaled <- multiplier * scaled$z
    validation <- bpe_validate_design(
      y ~ x - 1 | z_scaled - 1, scaled, design,
      bpe_min_n_S = 500, bpe_equiv_margin = 0.08,
      scale_instrument = "residual_sd"
    )
    fit <- spliv(
      y ~ x - 1 | z_scaled - 1, scaled, method = "bpe",
      bpe_design = design, bpe_min_n_S = 500,
      bpe_equiv_margin = 0.08, scale_instrument = "residual_sd"
    )
    list(validation = validation, fit = fit)
  }

  fits <- lapply(c(1, 100, 0.01), run_scaled)
  reference <- fits[[1]]
  for (candidate in fits[-1]) {
    expect_equal(
      candidate$validation$standardized_first_stage_effect,
      reference$validation$standardized_first_stage_effect,
      tolerance = 1e-9
    )
    expect_equal(
      candidate$validation$standardized_first_stage_ci,
      reference$validation$standardized_first_stage_ci,
      tolerance = 1e-9
    )
    expect_identical(candidate$validation$equivalence_passed,
                     reference$validation$equivalence_passed)
    expect_identical(candidate$validation$eligibility_passed,
                     reference$validation$eligibility_passed)
    expect_equal(candidate$fit$estimates, reference$fit$estimates, tolerance = 1e-8)
  }
  expect_true(reference$validation$eligibility_passed)
  expect_equal(reference$validation$raw_first_stage_coefficient,
               reference$validation$first_stage_coefficient)
  expect_equal(reference$validation$raw_first_stage_ci,
               reference$validation$first_stage_ci)
  expect_identical(
    reference$validation$equivalence_scale,
    "residual_treatment_sd_per_residual_instrument_sd"
  )
})

test_that("raw BPE equivalence scale compares the raw first-stage interval", {
  d <- make_bpe_redesign_data(n = 5000, pi_S = 0.02, pi_notS = 1, seed = 3104)
  design <- bpe_design(
    "Inactive subset", ~ inactive_region,
    rationale = "The treatment channel is absent in this pre-specified subset.",
    transportability_rationale = "The direct-effect mechanism is assumed to apply to the target sample."
  )
  validation <- bpe_validate_design(
    y ~ x - 1 | z - 1, d, design,
    bpe_min_n_S = 500, bpe_equiv_margin = 0.08,
    scale_instrument = "none"
  )
  expect_identical(validation$equivalence_scale, "raw_first_stage_coefficient")
  expect_equal(validation$equivalence_ci, validation$raw_first_stage_ci)
  ci <- validation$raw_first_stage_ci["z", ]
  expect_identical(
    validation$equivalence_passed,
    isTRUE(ci[["lower"]] >= -0.08 && ci[["upper"]] <= 0.08)
  )
})

test_that("transportability rationale is a confirmatory eligibility requirement", {
  d <- make_bpe_redesign_data(n = 5000, pi_S = 0, pi_notS = 1, seed = 3105)
  make_design <- function(value) {
    bpe_design(
      "Inactive subset", ~ inactive_region,
      rationale = "The treatment channel is absent in this pre-specified subset.",
      transportability_rationale = value
    )
  }

  for (bad_value in list(NULL, "", "   ")) {
    design <- make_design(bad_value)
    validation <- bpe_validate_design(
      y ~ x - 1 | z - 1, d, design,
      bpe_min_n_S = 500, bpe_equiv_margin = 0.08
    )
    expect_false(validation$transportability_rationale_passed)
    expect_false(validation$eligibility_passed)
    expect_match(validation$message, "transportability_rationale", fixed = TRUE)

    fit <- spliv(
      y ~ x - 1 | z - 1, d, method = "bpe",
      bpe_design = design, bpe_min_n_S = 500,
      bpe_equiv_margin = 0.08, bpe_not_applicable = "na"
    )
    expect_true(all(is.na(fit$estimates$estimate)))
  }

  valid <- bpe_validate_design(
    y ~ x - 1 | z - 1, d,
    make_design("The direct-effect mechanism is assumed to apply to the target sample."),
    bpe_min_n_S = 500, bpe_equiv_margin = 0.08
  )
  expect_true(valid$transportability_rationale_passed)
  expect_true(valid$eligibility_passed)
})

test_that("BPE transport modes and standalone validation match final estimation", {
  d <- make_bpe_redesign_data(n = 6000, pi_S = 0, pi_notS = 1, seed = 3106)
  design <- bpe_design(
    "Inactive subset", ~ inactive_region,
    rationale = "The treatment channel is absent in this pre-specified subset.",
    transportability_rationale = "The direct-effect mechanism is assumed to apply to the target sample."
  )
  sampling <- bpe_validate_design(
    y ~ x - 1 | z - 1, d, design,
    bpe_min_n_S = 500, bpe_equiv_margin = 0.08
  )
  conservative <- bpe_validate_design(
    y ~ x - 1 | z - 1, d, design,
    bpe_min_n_S = 500, bpe_equiv_margin = 0.08,
    bpe_transport = "conservative", bpe_transport_kappa = 0.5
  )
  fit <- spliv(
    y ~ x - 1 | z - 1, d, method = "bpe",
    bpe_design = design, bpe_min_n_S = 500,
    bpe_equiv_margin = 0.08
  )

  expect_identical(sampling$transport_mode, "sampling")
  expect_identical(conservative$transport_mode, "conservative")
  expect_equal(conservative$transport_covariance,
               1.5 * sampling$transport_covariance, tolerance = 1e-12)
  expect_equal(conservative$prior_mu_full, sampling$prior_mu_full, tolerance = 1e-12)
  expect_equal(fit$mu_used, unname(sampling$prior_mu_full), tolerance = 1e-12)
  expect_equal(fit$Omega_used, sampling$prior_Omega_full, tolerance = 1e-12)
})

test_that("BPE reports insufficient subset information without estimating", {
  d <- make_bpe_redesign_data(n = 400, pi_S = 0, pi_notS = 1, seed = 3107)
  design <- bpe_design(
    "Inactive subset", ~ inactive_region,
    rationale = "The treatment channel is absent in this pre-specified subset.",
    transportability_rationale = "The direct-effect mechanism is assumed to apply to the target sample."
  )
  too_small <- bpe_validate_design(
    y ~ x - 1 | z - 1, d, design,
    bpe_min_n_S = 300, bpe_equiv_margin = 0.2
  )
  expect_false(too_small$eligibility_checks$minimum_n)
  expect_false(too_small$eligibility_passed)

  no_variation <- d
  no_variation$z[no_variation$inactive_region] <- 1
  variation <- bpe_validate_design(
    y ~ x - 1 | z - 1, no_variation, design,
    bpe_min_n_S = 100, bpe_equiv_margin = 0.2
  )
  expect_false(variation$eligibility_passed)
  expect_match(variation$message, "(?i)variation|rank|constant")
})

test_that("missingness, cluster inputs, weak instruments, and singular designs are handled clearly", {
  d <- make_synth_panel(n_gid = 20, n_t = 15, seed = 3108)
  d$y[c(1, 7)] <- NA_real_
  d$w1[11] <- NA_real_
  d$gid[19] <- NA
  keep <- complete.cases(d[, c("y", "x", "z", "w1", "w2", "gid")])

  formula_cluster <- spliv(
    y ~ x + w1 + w2 | z + w1 + w2, d,
    method = "uci", delta = 0, vcov = "cluster", cluster = ~ gid
  )
  clean <- make_synth_panel(n_gid = 20, n_t = 15, seed = 3108)
  character_cluster <- spliv(
    y ~ x + w1 + w2 | z + w1 + w2, clean,
    method = "uci", delta = 0, vcov = "cluster",
    cluster = as.character(clean$gid)
  )
  clean_formula_cluster <- spliv(
    y ~ x + w1 + w2 | z + w1 + w2, clean,
    method = "uci", delta = 0, vcov = "cluster", cluster = ~ gid
  )
  expect_equal(length(formula_cluster$internals$y), sum(keep))
  expect_equal(clean_formula_cluster$estimates, character_cluster$estimates, tolerance = 1e-12)

  set.seed(3109)
  weak <- data.frame(z = rnorm(500))
  weak$x <- 1e-5 * weak$z + rnorm(500)
  weak$y <- weak$x + rnorm(500)
  weak_fit <- spliv(y ~ x | z, weak, method = "uci", delta = 0)
  expect_true(all(is.finite(weak_fit$estimates$conf.low)))

  singular <- data.frame(y = rnorm(100), x = rnorm(100), z = 1, w = rnorm(100))
  expect_error(
    spliv(y ~ x + w | z + w, singular, method = "uci", delta = 0),
    "(?i)singular|rank|variation|residualized instrument SD"
  )
  singular$w2 <- singular$w
  singular$z <- rnorm(100)
  expect_error(
    spliv(y ~ x + w + w2 | z + w + w2, singular, method = "uci", delta = 0),
    "(?i)singular|rank"
  )
})

test_that("all plotting interfaces validate terms and return inputs invisibly", {
  d <- make_pattern_sensitivity_data(n = 300, seed = 3110)
  fit <- spliv(y ~ x + w | z + w, d, method = "uci", delta = 0)
  path <- spliv_sensitivity_path(
    y ~ x + w | z + w, d, method = "uci", delta_grid = c(0, 0.05)
  )
  low_level <- sp_sensitivity_uci_support(
    y ~ x + w | z + w, d, term = "x", inst_vary = "z",
    delta_grid = c(0, 0.05), grid = 3
  )

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({
    if (grDevices::dev.cur() > 1) grDevices::dev.off()
    unlink(plot_file)
  }, add = TRUE)

  expect_identical(plot(fit), fit)
  expect_identical(plot(path, term = "x"), path)
  expect_identical(plot_sp_sensitivity(low_level), low_level)
  expect_error(plot(path, term = "absent"), "not found")
})
