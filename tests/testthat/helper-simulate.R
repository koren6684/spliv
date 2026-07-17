make_synth_panel <- function(n_gid = 80, n_t = 40, seed = 1L) {
  set.seed(seed)
  n <- n_gid * n_t
  gid <- rep(seq_len(n_gid), each = n_t)
  t <- rep(seq_len(n_t), times = n_gid)
  month <- ((t - 1L) %% 12L) + 1L
  year <- 2000L + ((t - 1L) %/% 12L)
  ym <- paste0(year, "-", sprintf("%02d", month))

  g_fe <- rnorm(n_gid, sd = 1.0)[gid]
  ym_levels <- unique(ym)
  ym_id <- match(ym, ym_levels)
  ym_fe <- rnorm(length(ym_levels), sd = 0.4)[ym_id]

  w1 <- rnorm(n)
  w2 <- rnorm(n)
  z <- 0.8 * w1 + rnorm(n)
  u <- rnorm(n)
  x <- 0.7 * z + 0.4 * w2 + g_fe + ym_fe + u
  eps <- rnorm(n)
  y <- 1.2 * x + 0.3 * w1 - 0.2 * w2 + 0.1 * z + g_fe + ym_fe + eps

  data.frame(
    y = y,
    x = x,
    z = z,
    w1 = w1,
    w2 = w2,
    gid = factor(gid),
    month = factor(month),
    year = factor(year),
    ym = factor(ym)
  )
}

make_bpe_redesign_data <- function(n = 8000,
                                   pi_S = 0,
                                   pi_notS = 1.1,
                                   beta = 1.2,
                                   gamma = 0.4,
                                   seed = 1L) {
  set.seed(seed)
  z <- rnorm(n)
  u <- rnorm(n)
  e <- rnorm(n)
  inactive_region <- rep(c(TRUE, FALSE), each = n / 2)
  pi_vec <- ifelse(inactive_region, pi_S, pi_notS)
  x <- pi_vec * z + u
  y <- beta * x + gamma * z + e
  data.frame(
    y = y,
    x = x,
    z = z,
    inactive_region = inactive_region,
    theoretical_condition = as.integer(inactive_region),
    treatment_channel_absent = inactive_region,
    logical_subset = inactive_region
  )
}

make_pattern_sensitivity_data <- function(n = 6000,
                                          beta = 1.2,
                                          pi = 1.0,
                                          gamma = 0.3,
                                          seed = 1L) {
  set.seed(seed)
  z <- rnorm(n)
  w <- rnorm(n)
  treatment_channel_exposure <- stats::pnorm(w)
  high_exposure_region <- treatment_channel_exposure > stats::median(treatment_channel_exposure)
  x <- pi * z + 0.5 * w + rnorm(n)
  y <- beta * x + gamma * treatment_channel_exposure * z + 0.4 * w + rnorm(n)

  data.frame(
    y = y,
    x = x,
    z = z,
    w = w,
    treatment_channel_exposure = treatment_channel_exposure,
    high_exposure_region = high_exposure_region,
    logical_pattern = high_exposure_region
  )
}
