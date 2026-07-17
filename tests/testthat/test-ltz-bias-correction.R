test_that("LTZ scalar correction matches mu/pi in just-identified model", {
  set.seed(20260220)

  n <- 40000
  beta_true <- 1.25
  pi_true <- 1.5
  mu_true <- 0.6

  z <- rnorm(n)
  u <- rnorm(n)
  e <- rnorm(n)

  x <- pi_true * z + u
  y <- beta_true * x + mu_true * z + e

  X <- matrix(x, ncol = 1, dimnames = list(NULL, "x"))
  Z <- matrix(z, ncol = 1, dimnames = list(NULL, "z"))

  iv <- spliv:::.iv_2sls_mats(y = y, X = X, Z = Z, vcov = "iid")
  ltz <- spliv:::conley_ltz_mats(
    y = y,
    X = X,
    Z = Z,
    mu = mu_true,
    Omega = matrix(0, 1, 1),
    vcov = "iid",
    coef_names = "x",
    inst_names = "z"
  )

  expected_shift <- mu_true / pi_true
  observed_shift <- as.numeric(iv$beta - ltz$beta)
  naive_bias <- as.numeric(iv$beta - beta_true)
  corrected_bias <- as.numeric(ltz$beta - beta_true)

  expect_equal(observed_shift, expected_shift, tolerance = 0.05)
  expect_equal(naive_bias, expected_shift, tolerance = 0.05)
  expect_lt(abs(corrected_bias), 0.05)
})

test_that("LTZ supports a generalized direct-effect regressor for patterned sensitivity", {
  set.seed(20260424)

  n <- 40000
  beta_true <- 1.1
  pi_true <- 1.4
  theta_true <- 0.5

  z <- rnorm(n)
  exposure <- stats::pnorm(rnorm(n))
  g <- exposure * z
  u <- rnorm(n)
  e <- rnorm(n)

  x <- pi_true * z + u
  y <- beta_true * x + theta_true * g + e

  X <- matrix(x, ncol = 1, dimnames = list(NULL, "x"))
  Z <- matrix(z, ncol = 1, dimnames = list(NULL, "z"))
  G <- matrix(g, ncol = 1, dimnames = list(NULL, "theta_z_pattern"))

  ltz_pattern <- spliv:::sp_ltz_mats(
    y = y,
    X = X,
    Z = Z,
    mu = theta_true,
    Omega = matrix(0, 1, 1),
    vcov = "iid",
    coef_names = "x",
    inst_names = "z",
    direct_effect = G,
    direct_effect_names = "theta_z_pattern"
  )

  iv_adjusted <- spliv:::.iv_2sls_mats(y = y - theta_true * g, X = X, Z = Z, vcov = "iid")
  expect_equal(as.numeric(ltz_pattern$beta), as.numeric(iv_adjusted$beta), tolerance = 1e-10)
})
