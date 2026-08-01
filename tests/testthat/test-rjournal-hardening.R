test_that("default spliv() is conventional UCI at delta zero", {
  set.seed(3001)
  n <- 500
  z <- rnorm(n)
  x <- 0.8 * z + rnorm(n)
  y <- 1.3 * x + rnorm(n)
  d <- data.frame(y, x, z)

  default_fit <- spliv(y ~ x | z, d, vcov = "iid")
  explicit_fit <- spliv(y ~ x | z, d, method = "uci", delta = 0, vcov = "iid")

  expect_identical(default_fit$method, "uci")
  expect_equal(default_fit$estimates, explicit_fit$estimates, tolerance = 1e-12)
  expect_equal(default_fit$grid$delta, 0)
})

test_that("UCI delta-zero agrees with fixest IV with and without fixed effects", {
  set.seed(3002)
  n <- 500
  z <- rnorm(n)
  x <- 0.8 * z + rnorm(n)
  y <- 1.3 * x + rnorm(n)
  d <- data.frame(y, x, z)
  fit <- spliv(y ~ x | z, d, method = "uci", delta = 0, vcov = "iid")
  ref <- fixest::feols(y ~ 1 | 0 | x ~ z, data = d, vcov = "iid")
  fit_beta <- with(fit$estimates[fit$estimates$term == "x", ], (conf.low + conf.high) / 2)
  expect_equal(as.numeric(fit_beta), unname(stats::coef(ref)[["fit_x"]]), tolerance = 1e-10)

  panel <- make_synth_panel(n_gid = 12, n_t = 24, seed = 3003)
  f <- y ~ x + w1 + w2 | z + w1 + w2
  fit_fe <- spliv(f, panel, fe = ~ gid + ym, fe_engine = "fixest",
                  method = "uci", delta = 0, vcov = "iid")
  ref_fe <- fixest::feols(y ~ w1 + w2 | gid + ym | x ~ z, data = panel, vcov = "iid")
  fit_beta_fe <- with(fit_fe$estimates[fit_fe$estimates$term == "x", ],
                      (conf.low + conf.high) / 2)
  expect_equal(as.numeric(fit_beta_fe), unname(stats::coef(ref_fe)[["fit_x"]]),
               tolerance = 1e-10)
})

test_that("nonzero UCI adjustment matches an adjusted-outcome IV point estimate", {
  set.seed(3004)
  n <- 500
  z <- rnorm(n)
  x <- 0.8 * z + rnorm(n)
  mu <- 0.2
  y <- 1.3 * x + mu * z + rnorm(n)
  d <- data.frame(y, x, z)
  fit <- spliv(y ~ x | z, d, method = "uci", scale_instrument = "none",
               grid = list(gmin = mu, gmax = mu, steps = 1), vcov = "iid")
  d$adjusted_y <- d$y - mu * d$z
  ref <- fixest::feols(adjusted_y ~ 1 | 0 | x ~ z, data = d, vcov = "iid")
  fit_beta <- with(fit$estimates[fit$estimates$term == "x", ], (conf.low + conf.high) / 2)
  expect_equal(as.numeric(fit_beta), unname(stats::coef(ref)[["fit_x"]]), tolerance = 1e-10)
  # Confidence widths can differ at finite samples because spliv uses a
  # normal critical value while fixest::confint() uses its own small-sample
  # convention; the point estimate is the numerical equivalence target.
})

test_that("lfe and fixest fixed-effect backends agree when lfe is installed", {
  testthat::skip_if_not_installed("lfe")
  panel <- make_synth_panel(n_gid = 10, n_t = 24, seed = 3005)
  f <- y ~ x + w1 + w2 | z + w1 + w2
  fixest_fit <- spliv(f, panel, fe = ~ gid + ym, fe_engine = "fixest",
                      method = "uci", delta = 0, vcov = "iid")
  lfe_fit <- spliv(f, panel, fe = ~ gid + ym, fe_engine = "lfe",
                   method = "uci", delta = 0, vcov = "iid")
  expect_equal(fixest_fit$estimates, lfe_fit$estimates, tolerance = 1e-10)
})

test_that("exploratory subset diagnostics report F but do not select on it", {
  set.seed(3006)
  d <- data.frame(y = rnorm(80), x = rnorm(80), z = rnorm(80),
                  inactive = rep(c(TRUE, FALSE), each = 40))
  out <- NULL
  expect_warning(
    out <- bpe_explore_subsets(d, y ~ x | z, rules = list(inactive = ~ inactive)),
    "exploratory"
  )
  expect_true("F_S" %in% names(out))
  expect_false("screen_F_ok" %in% names(out))
  expect_true(is.finite(out$F_S[[1]]) || is.na(out$F_S[[1]]))
})

test_that("low-level sensitivity helpers and plotting validate terms", {
  set.seed(3007)
  d <- data.frame(y = rnorm(80), x = rnorm(80), z = rnorm(80))
  out <- sp_sensitivity_ltz_uniform01_as_normal(
    y ~ x | z, d, term = "x", inst_vary = "z",
    delta_grid = c(0, 0.1), scale_instrument = "none"
  )
  expect_true(is.data.frame(out))
  expect_equal(unique(out$method), "LTZ (Normal approx to U(0,delta))")

  path <- spliv_sensitivity_path(y ~ x | z, d, method = "uci",
                                 delta_grid = c(0, 0.1), vcov = "iid")
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({
    if (grDevices::dev.cur() > 1) grDevices::dev.off()
    unlink(plot_file)
  }, add = TRUE)
  expect_silent(plot_sp_sensitivity(path, term = "x"))
  expect_error(plot_sp_sensitivity(path, term = "missing"), "not found")
})
