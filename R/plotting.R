.resolve_plot_inputs <- function(formula_or_fit, data = NULL) {
  if (inherits(formula_or_fit, "plausexog_fit")) {
    fit <- formula_or_fit
    return(list(
      y = fit$internals$y,
      X = fit$internals$X,
      Z = fit$internals$Z,
      W = fit$internals$W,
      cluster_id = fit$internals$cluster_id,
      vcov = fit$vcov,
      scale_instrument = fit$scale_instrument %||% "residual_sd",
      residualized_instrument_sd = fit$residualized_instrument_sd %||%
        .instrument_residual_sds(
          X = fit$internals$X,
          Z = fit$internals$Z,
          W = fit$internals$W,
          fe_present = FALSE
        ),
      coef_names = colnames(fit$internals$X),
      inst_names = colnames(fit$internals$Z)
    ))
  }

  if (is.null(data)) {
    stop("`data` is required when first argument is a formula.")
  }
  parsed <- .iv_parse(formula_or_fit, data)
  W <- .build_W_from_XZ(parsed$X, parsed$Z)
  list(
    y = parsed$y,
    X = parsed$X,
    Z = parsed$Z,
    W = W,
    cluster_id = NULL,
    vcov = "hc1",
    scale_instrument = "residual_sd",
    residualized_instrument_sd = .instrument_residual_sds(
      X = parsed$X,
      Z = parsed$Z,
      W = W,
      fe_present = FALSE
    ),
    coef_names = parsed$coef_names,
    inst_names = parsed$inst_names
  )
}

.plot_scale_delta <- function(values, inst_vary, residualized_instrument_sd, scale_instrument) {
  values <- as.numeric(values)
  if (identical(scale_instrument, "none")) {
    return(values)
  }
  sd_inst <- residualized_instrument_sd[inst_vary]
  if (any(!is.finite(sd_inst)) || any(sd_inst <= 0)) {
    stop(
      "Residualized instrument SD must be positive when `scale_instrument = 'residual_sd'`."
    )
  }
  values / sd_inst
}

#' Build LTZ Prior Matrices for Chosen Instruments
#'
#' @param formula Formula or `plausexog_fit`.
#' @param data Data frame when `formula` is a formula.
#' @param inst_vary Instrument(s) with non-degenerate prior.
#' @param mean Prior mean(s).
#' @param sd Prior sd(s).
#' @param vcov Vcov type.
#' @param cluster Optional cluster ids.
#'
#' @return List with `mu`, `omega`, and instrument names.
#' @examples
#' set.seed(9)
#' d <- data.frame(y = rnorm(60), x = rnorm(60), z = rnorm(60))
#' sp_prior_ltz(y ~ x | z, d, inst_vary = "z", sd = 0.1)
#' @export
sp_prior_ltz <- function(formula, data = NULL,
                             inst_vary,
                             mean = 0,
                             sd = 1,
                             vcov = c("hc1", "hc0", "iid", "cluster"),
                             cluster = NULL) {
  vcov <- match.arg(tolower(vcov), c("hc1", "hc0", "iid", "cluster"))
  inp <- .resolve_plot_inputs(formula, data)

  inst_vary <- as.character(inst_vary)
  if (!all(inst_vary %in% inp$inst_names)) {
    stop("inst_vary not found in instrument matrix Z.")
  }
  idx <- match(inst_vary, inp$inst_names)

  if (length(mean) == 1) mean <- rep(mean, length(idx))
  if (length(sd) == 1) sd <- rep(sd, length(idx))

  mu <- rep(0, length(inp$inst_names))
  mu[idx] <- mean
  omega <- matrix(0, length(inp$inst_names), length(inp$inst_names))
  omega[idx, idx] <- diag(sd^2, nrow = length(idx))

  list(mu = mu, omega = omega, inst_names = inp$inst_names)
}

#' @rdname sp_prior_ltz
#' @export
conley_prior_ltz <- function(formula, data = NULL,
                             inst_vary,
                             mean = 0,
                             sd = 1,
                             vcov = c("hc1", "hc0", "iid", "cluster"),
                             cluster = NULL) {
  sp_prior_ltz(
    formula = formula,
    data = data,
    inst_vary = inst_vary,
    mean = mean,
    sd = sd,
    vcov = vcov,
    cluster = cluster
  )
}

#' LTZ Sensitivity over Delta Grid
#'
#' @param formula Formula or `plausexog_fit`.
#' @param data Data frame when `formula` is a formula.
#' @param term Coefficient name to track.
#' @param inst_vary Instrument(s) with plausible violation.
#' @param delta_grid Delta grid.
#' @param mean_fun Function mapping delta to prior mean.
#' @param sd_fun Function mapping delta to prior sd.
#' @param level Confidence level.
#' @param vcov Vcov type.
#' @param cluster Optional cluster ids.
#' @param scale_instrument One of `"residual_sd"` (default) or `"none"`.
#'
#' @return Data frame with sensitivity path.
#' @export
sp_sensitivity_ltz_normal <- function(formula, data = NULL,
                                          term,
                                          inst_vary,
                                          delta_grid,
                                          mean_fun = function(delta) 0,
                                          sd_fun = function(delta) delta,
                                          level = 0.95,
                                          vcov = c("hc1", "hc0", "iid", "cluster"),
                                          cluster = NULL,
                                          scale_instrument = c("residual_sd", "none")) {
  vcov <- match.arg(tolower(vcov), c("hc1", "hc0", "iid", "cluster"))
  scale_instrument <- match.arg(scale_instrument)
  inp <- .resolve_plot_inputs(formula, data)
  if (inherits(formula, "plausexog_fit")) {
    vcov <- inp$vcov
    if (missing(scale_instrument)) {
      scale_instrument <- inp$scale_instrument
    }
  }

  out <- vector("list", length(delta_grid))
  for (i in seq_along(delta_grid)) {
    d <- delta_grid[i]
    mean_val <- .plot_scale_delta(
      values = mean_fun(d),
      inst_vary = inst_vary,
      residualized_instrument_sd = inp$residualized_instrument_sd,
      scale_instrument = scale_instrument
    )
    sd_val <- .plot_scale_delta(
      values = sd_fun(d),
      inst_vary = inst_vary,
      residualized_instrument_sd = inp$residualized_instrument_sd,
      scale_instrument = scale_instrument
    )
    prior <- sp_prior_ltz(
      formula = formula,
      data = data,
      inst_vary = inst_vary,
      mean = mean_val,
      sd = sd_val,
      vcov = vcov,
      cluster = cluster
    )

    ltz <- sp_ltz_mats(
      y = inp$y,
      X = inp$X,
      Z = inp$Z,
      mu = prior$mu,
      Omega = prior$omega,
      level = level,
      vcov = vcov,
      cluster_id = inp$cluster_id,
      coef_names = inp$coef_names,
      inst_names = inp$inst_names
    )$table

    row <- ltz[ltz$term == term, , drop = FALSE]
    if (nrow(row) != 1) stop("term not found or duplicated: ", term)

    out[[i]] <- data.frame(
      delta = d,
      estimate = row$estimate,
      conf.low = row$conf.low,
      conf.high = row$conf.high,
      method = "LTZ (Normal prior)",
      row.names = NULL
    )
  }
  do.call(rbind, out)
}

#' @rdname sp_sensitivity_ltz_normal
#' @export
conley_sensitivity_ltz_normal <- function(formula, data = NULL,
                                          term,
                                          inst_vary,
                                          delta_grid,
                                          mean_fun = function(delta) 0,
                                          sd_fun = function(delta) delta,
                                          level = 0.95,
                                          vcov = c("hc1", "hc0", "iid", "cluster"),
                                          cluster = NULL,
                                          scale_instrument = c("residual_sd", "none")) {
  sp_sensitivity_ltz_normal(
    formula = formula,
    data = data,
    term = term,
    inst_vary = inst_vary,
    delta_grid = delta_grid,
    mean_fun = mean_fun,
    sd_fun = sd_fun,
    level = level,
    vcov = vcov,
    cluster = cluster,
    scale_instrument = scale_instrument
  )
}

#' UCI Sensitivity over Delta Grid
#'
#' @param formula Formula or `plausexog_fit`.
#' @param data Data frame when `formula` is a formula.
#' @param term Coefficient name to track.
#' @param inst_vary Instrument(s) whose gamma bounds vary.
#' @param delta_grid Delta grid.
#' @param gmin_fun Function mapping delta to lower bound.
#' @param gmax_fun Function mapping delta to upper bound.
#' @param grid Points per instrument bound.
#' @param level Confidence level.
#' @param vcov Vcov type.
#' @param cluster Optional cluster ids.
#' @param scale_instrument One of `"residual_sd"` (default) or `"none"`.
#'
#' @return Data frame with sensitivity path.
#' @export
sp_sensitivity_uci_support <- function(formula, data = NULL,
                                           term,
                                           inst_vary,
                                           delta_grid,
                                           gmin_fun = function(delta) -delta,
                                           gmax_fun = function(delta) delta,
                                           grid = 41,
                                           level = 0.95,
                                           vcov = c("hc1", "hc0", "iid", "cluster"),
                                           cluster = NULL,
                                           scale_instrument = c("residual_sd", "none")) {
  vcov <- match.arg(tolower(vcov), c("hc1", "hc0", "iid", "cluster"))
  scale_instrument <- match.arg(scale_instrument)
  inp <- .resolve_plot_inputs(formula, data)
  if (inherits(formula, "plausexog_fit")) {
    vcov <- inp$vcov
    if (missing(scale_instrument)) {
      scale_instrument <- inp$scale_instrument
    }
  }

  inst_idx <- match(inst_vary, inp$inst_names)
  if (any(is.na(inst_idx))) {
    stop("Some `inst_vary` names were not found in Z.")
  }

  out <- vector("list", length(delta_grid))
  for (i in seq_along(delta_grid)) {
    d <- delta_grid[i]
    p <- length(inst_idx)
    gmin <- .plot_scale_delta(
      values = gmin_fun(d),
      inst_vary = inst_vary,
      residualized_instrument_sd = inp$residualized_instrument_sd,
      scale_instrument = scale_instrument
    )
    gmax <- .plot_scale_delta(
      values = gmax_fun(d),
      inst_vary = inst_vary,
      residualized_instrument_sd = inp$residualized_instrument_sd,
      scale_instrument = scale_instrument
    )
    if (length(gmin) == 1) gmin <- rep(gmin, p)
    if (length(gmax) == 1) gmax <- rep(gmax, p)
    steps <- if (length(grid) == 1) rep(grid, p) else grid

    gamma_grid <- .make_gamma_grid(gmin = gmin, gmax = gmax, steps = steps)
    ci <- sp_uci_mats(
      y = inp$y,
      X = inp$X,
      Z = inp$Z,
      inst_idx = inst_idx,
      gamma_grid = gamma_grid,
      level = level,
      vcov = vcov,
      cluster_id = inp$cluster_id,
      coef_names = inp$coef_names
    )

    row <- ci[ci$term == term, , drop = FALSE]
    if (nrow(row) != 1) stop("term not found or duplicated: ", term)

    out[[i]] <- data.frame(
      delta = d,
      estimate = NA_real_,
      conf.low = row$conf.low,
      conf.high = row$conf.high,
      method = "UCI (Support only)",
      row.names = NULL
    )
  }

  do.call(rbind, out)
}

#' @rdname sp_sensitivity_uci_support
#' @export
conley_sensitivity_uci_support <- function(formula, data = NULL,
                                           term,
                                           inst_vary,
                                           delta_grid,
                                           gmin_fun = function(delta) -delta,
                                           gmax_fun = function(delta) delta,
                                           grid = 41,
                                           level = 0.95,
                                           vcov = c("hc1", "hc0", "iid", "cluster"),
                                           cluster = NULL,
                                           scale_instrument = c("residual_sd", "none")) {
  sp_sensitivity_uci_support(
    formula = formula,
    data = data,
    term = term,
    inst_vary = inst_vary,
    delta_grid = delta_grid,
    gmin_fun = gmin_fun,
    gmax_fun = gmax_fun,
    grid = grid,
    level = level,
    vcov = vcov,
    cluster = cluster,
    scale_instrument = scale_instrument
  )
}

#' LTZ Sensitivity with Normal Approximation to U(0, delta)
#'
#' @inheritParams sp_sensitivity_ltz_normal
#' @return Data frame with sensitivity path.
#' @export
sp_sensitivity_ltz_uniform01_as_normal <- function(formula, data = NULL,
                                                       term, inst_vary, delta_grid,
                                                       level = 0.95,
                                                       vcov = c("hc1", "hc0", "iid", "cluster"),
                                                       cluster = NULL,
                                                       scale_instrument = c("residual_sd", "none")) {
  out <- sp_sensitivity_ltz_normal(
    formula = formula,
    data = data,
    term = term,
    inst_vary = inst_vary,
    delta_grid = delta_grid,
    mean_fun = function(d) d / 2,
    sd_fun = function(d) sqrt(d^2 / 12),
    level = level,
    vcov = vcov,
    cluster = cluster,
    scale_instrument = scale_instrument
  )
  out$method <- "LTZ (Normal approx to U(0,delta))"
  out
}

#' @rdname sp_sensitivity_ltz_uniform01_as_normal
#' @export
conley_sensitivity_ltz_uniform01_as_normal <- function(formula, data = NULL,
                                                       term, inst_vary, delta_grid,
                                                       level = 0.95,
                                                       vcov = c("hc1", "hc0", "iid", "cluster"),
                                                       cluster = NULL,
                                                       scale_instrument = c("residual_sd", "none")) {
  sp_sensitivity_ltz_uniform01_as_normal(
    formula = formula,
    data = data,
    term = term,
    inst_vary = inst_vary,
    delta_grid = delta_grid,
    level = level,
    vcov = vcov,
    cluster = cluster,
    scale_instrument = scale_instrument
  )
}

.validate_sensitivity_delta_grid <- function(delta_grid) {
  delta_grid <- as.numeric(delta_grid)
  if (!length(delta_grid)) {
    stop("`delta_grid` must contain at least one value.")
  }
  if (any(!is.finite(delta_grid)) || any(delta_grid < 0)) {
    stop("`delta_grid` must contain non-negative finite values.")
  }
  delta_grid
}

.check_sensitivity_path_dots <- function(dots) {
  reserved <- c("formula", "data", "method", "delta", "violation_pattern")
  bad_reserved <- intersect(names(dots), reserved)
  if (length(bad_reserved)) {
    stop(
      "`spliv_sensitivity_path()` controls ",
      paste(sprintf("`%s`", bad_reserved), collapse = ", "),
      "; remove those arguments from `...`."
    )
  }

  if (!is.null(dots$prior)) {
    stop("`spliv_sensitivity_path()` varies `delta` directly; do not supply `prior` in `...`.")
  }
  if (isTRUE(dots$bpe %||% FALSE)) {
    stop("`spliv_sensitivity_path()` currently supports LTZ and UCI only; BPE is intentionally excluded.")
  }
  if (!is.null(dots$grid)) {
    if (!is.list(dots$grid)) {
      stop("`grid` must be a list when supplied through `...`.")
    }
    bad_grid <- intersect(names(dots$grid), c("delta", "gmin", "gmax"))
    if (length(bad_grid)) {
      stop(
        "`spliv_sensitivity_path()` controls the delta/theta range; remove ",
        paste(sprintf("`grid$%s`", bad_grid), collapse = ", "),
        " from `...`."
      )
    }
  }
}

.fit_for_sensitivity_delta <- function(formula, data, method, delta, violation_pattern, dots) {
  do.call(
    spliv,
    c(
      list(
        formula = formula,
        data = data,
        method = method,
        delta = delta,
        violation_pattern = violation_pattern
      ),
      dots
    )
  )
}

.extract_interval_column <- function(estimates, names_to_try) {
  for (nm in names_to_try) {
    if (nm %in% names(estimates)) {
      return(as.numeric(estimates[[nm]]))
    }
  }
  NULL
}

.baseline_reference_rows <- function(fit, method) {
  .extract_sensitivity_rows(
    fit = fit,
    delta = 0,
    method = method,
    baseline_ref = NULL,
    error_message = NA_character_
  )
}

.extract_sensitivity_rows <- function(fit,
                                      delta,
                                      method,
                                      baseline_ref = NULL,
                                      error_message = NA_character_) {
  estimates <- fit$estimates
  if (!is.data.frame(estimates)) {
    stop("`fit$estimates` must be a data frame.")
  }

  term <- estimates$term %||% colnames(fit$internals$X)
  conf_low <- .extract_interval_column(estimates, c("conf.low", "conf_low"))
  conf_high <- .extract_interval_column(estimates, c("conf.high", "conf_high"))
  if (is.null(conf_low) || is.null(conf_high)) {
    stop("Sensitivity extraction requires confidence interval columns in `fit$estimates`.")
  }

  raw_estimate <- .extract_interval_column(estimates, c("estimate"))
  if (is.null(raw_estimate) && length(fit$beta_hat) == length(term)) {
    raw_estimate <- as.numeric(fit$beta_hat)
  }
  if (is.null(raw_estimate)) {
    raw_estimate <- rep(NA_real_, length(term))
  }

  raw_se <- .extract_interval_column(estimates, c("std.error", "std_error", "se"))
  if (is.null(raw_se)) {
    raw_se <- rep(NA_real_, length(term))
  }

  midpoint <- ifelse(
    is.finite(conf_low) & is.finite(conf_high),
    0.5 * (conf_low + conf_high),
    NA_real_
  )

  if (is.null(baseline_ref)) {
    baseline_estimate <- raw_estimate
    baseline_estimate[!is.finite(baseline_estimate)] <- midpoint[!is.finite(baseline_estimate)]
    baseline_conf_low <- conf_low
    baseline_conf_high <- conf_high
  } else {
    idx <- match(term, baseline_ref$term)
    baseline_estimate <- baseline_ref$baseline_estimate[idx]
    baseline_conf_low <- baseline_ref$baseline_conf_low[idx]
    baseline_conf_high <- baseline_ref$baseline_conf_high[idx]
  }

  estimate_out <- raw_estimate
  estimate_out[!is.finite(estimate_out)] <- baseline_estimate[!is.finite(estimate_out)]

  contains_zero <- ifelse(
    is.finite(conf_low) & is.finite(conf_high),
    conf_low <= 0 & conf_high >= 0,
    NA
  )
  significant_at_level <- ifelse(is.na(contains_zero), NA, !contains_zero)

  crosses_baseline_sign <- rep(NA, length(term))
  positive_baseline <- is.finite(baseline_estimate) & baseline_estimate > 0
  negative_baseline <- is.finite(baseline_estimate) & baseline_estimate < 0
  crosses_baseline_sign[positive_baseline] <- conf_low[positive_baseline] <= 0
  crosses_baseline_sign[negative_baseline] <- conf_high[negative_baseline] >= 0

  pattern_info <- fit$violation_pattern
  pattern_name <- if (is.null(pattern_info)) {
    "Uniform direct effect"
  } else {
    pattern_info$name
  }
  pattern_type <- if (is.null(pattern_info)) {
    "uniform"
  } else {
    pattern_info$pattern_type %||% NA_character_
  }

  data.frame(
    term = as.character(term),
    delta = as.numeric(delta),
    method = as.character(method),
    estimate = as.numeric(estimate_out),
    conf_low = as.numeric(conf_low),
    conf_high = as.numeric(conf_high),
    contains_zero = as.logical(contains_zero),
    pattern_name = rep(pattern_name, length(term)),
    pattern_type = rep(pattern_type, length(term)),
    violation_pattern_used = rep(!is.null(pattern_info), length(term)),
    scale_instrument = rep(fit$scale_instrument %||% NA_character_, length(term)),
    nobs = rep(length(fit$internals$y %||% numeric(0)), length(term)),
    se = as.numeric(raw_se),
    theta_min = if (identical(method, "uci")) rep(-delta, length(term)) else rep(NA_real_, length(term)),
    theta_max = if (identical(method, "uci")) rep(delta, length(term)) else rep(NA_real_, length(term)),
    baseline_estimate = as.numeric(baseline_estimate),
    baseline_conf_low = as.numeric(baseline_conf_low),
    baseline_conf_high = as.numeric(baseline_conf_high),
    crosses_baseline_sign = as.logical(crosses_baseline_sign),
    significant_at_level = as.logical(significant_at_level),
    error = rep(as.character(error_message), length(term)),
    stringsAsFactors = FALSE
  )
}

.error_sensitivity_rows <- function(baseline_ref, delta, method, error_message) {
  out <- baseline_ref
  out$delta <- as.numeric(delta)
  out$method <- as.character(method)
  out$estimate <- NA_real_
  out$conf_low <- NA_real_
  out$conf_high <- NA_real_
  out$contains_zero <- NA
  out$se <- NA_real_
  out$theta_min <- if (identical(method, "uci")) -delta else NA_real_
  out$theta_max <- if (identical(method, "uci")) delta else NA_real_
  out$crosses_baseline_sign <- NA
  out$significant_at_level <- NA
  out$error <- rep(as.character(error_message), nrow(out))
  out
}

.compute_tipping_points <- function(path) {
  terms <- unique(path$term)
  tipping <- stats::setNames(rep(NA_real_, length(terms)), terms)
  messages <- stats::setNames(rep("", length(terms)), terms)

  for (term_i in terms) {
    rows_i <- path[path$term == term_i, , drop = FALSE]
    rows_i <- rows_i[order(rows_i$delta), , drop = FALSE]

    baseline_contains_zero <- with(
      rows_i[1, , drop = FALSE],
      is.finite(baseline_conf_low) && is.finite(baseline_conf_high) &&
        baseline_conf_low <= 0 && baseline_conf_high >= 0
    )

    if (isTRUE(baseline_contains_zero)) {
      tipping[[term_i]] <- 0
      messages[[term_i]] <- "Baseline interval already includes zero."
      next
    }

    cross_idx <- which(isTRUE(rows_i$contains_zero) | rows_i$contains_zero %in% TRUE)
    if (!length(cross_idx)) {
      tipping[[term_i]] <- NA_real_
      messages[[term_i]] <- "No zero crossing occurred on the supplied delta grid."
      next
    }

    tipping[[term_i]] <- min(rows_i$delta[cross_idx], na.rm = TRUE)
  }

  list(values = tipping, messages = messages)
}

#' Sensitivity Path over Delta Grid
#'
#' Runs `spliv()` repeatedly over a user-supplied `delta_grid` and returns a
#' tidy sensitivity-path object. This is the recommended workflow for patterned
#' or uniform sensitivity analysis: users should usually report how intervals
#' change over a range of `delta` values rather than selecting one arbitrary
#' sensitivity level.
#'
#' @param formula IV formula `y ~ X | Z`.
#' @param data Data frame.
#' @param method One of `"uci"` or `"ltz"`. `spliv_sensitivity_path()`
#'   intentionally excludes confirmatory BPE.
#' @param delta_grid Non-negative sensitivity grid. For UCI, each value `d`
#'   implies theta bounds `[-d, +d]`.
#' @param violation_pattern Optional `spliv_pattern()` object. If omitted, the
#'   path uses the package's backward-compatible uniform direct-effect pattern.
#' @param stop_on_error Logical; if `TRUE` (default), stop on the first failed
#'   fit. If `FALSE`, record `NA` rows and store the error message in the
#'   returned `error` column.
#' @param ... Additional named arguments passed through to `spliv()`, such as
#'   `fe`, `vcov`, `cluster`, `scale_instrument`, or `grid = list(level = 0.95)`.
#'
#' @return A data frame with class `c("spliv_sensitivity_path", "data.frame")`.
#'   The returned object includes path columns such as `delta`, `method`,
#'   `estimate`, `conf_low`, `conf_high`, `contains_zero`, and pattern metadata,
#'   plus attributes containing the original call, the supplied grid, and a
#'   tipping-point summary.
#' @examples
#' set.seed(5)
#' d <- data.frame(y = rnorm(80), x = rnorm(80), z = rnorm(80))
#' p <- spliv_sensitivity_path(y ~ x | z, d, method = "uci",
#'   delta_grid = c(0, 0.1), vcov = "hc1")
#' head(p)
#' @export
spliv_sensitivity_path <- function(formula,
                                   data,
                                   method = c("uci", "ltz"),
                                   delta_grid = seq(0, 0.20, by = 0.01),
                                   violation_pattern = NULL,
                                   stop_on_error = TRUE,
                                   ...) {
  method <- match.arg(tolower(method), c("uci", "ltz", "bpe"))
  if (identical(method, "bpe")) {
    stop("`spliv_sensitivity_path()` currently supports LTZ and UCI only; BPE is intentionally excluded from sensitivity paths.")
  }
  if (!is.logical(stop_on_error) || length(stop_on_error) != 1 || is.na(stop_on_error)) {
    stop("`stop_on_error` must be TRUE or FALSE.")
  }

  delta_grid <- .validate_sensitivity_delta_grid(delta_grid)
  dots <- list(...)
  .check_sensitivity_path_dots(dots)

  baseline_fit <- tryCatch(
    .fit_for_sensitivity_delta(
      formula = formula,
      data = data,
      method = method,
      delta = 0,
      violation_pattern = violation_pattern,
      dots = dots
    ),
    error = function(e) e
  )
  if (inherits(baseline_fit, "error")) {
    stop("Failed to fit the baseline sensitivity model at delta = 0: ", conditionMessage(baseline_fit), call. = FALSE)
  }

  baseline_rows <- .baseline_reference_rows(baseline_fit, method = method)

  out_list <- vector("list", length(delta_grid))
  zero_tol <- sqrt(.Machine$double.eps)
  for (i in seq_along(delta_grid)) {
    delta_i <- delta_grid[[i]]
    fit_i <- if (abs(delta_i) < zero_tol) {
      baseline_fit
    } else {
      tryCatch(
        .fit_for_sensitivity_delta(
          formula = formula,
          data = data,
          method = method,
          delta = delta_i,
          violation_pattern = violation_pattern,
          dots = dots
        ),
        error = function(e) e
      )
    }

    if (inherits(fit_i, "error")) {
      if (isTRUE(stop_on_error)) {
        stop(
          "Sensitivity path fit failed at delta = ", format(delta_i),
          ": ", conditionMessage(fit_i),
          call. = FALSE
        )
      }
      out_list[[i]] <- .error_sensitivity_rows(
        baseline_ref = baseline_rows,
        delta = delta_i,
        method = method,
        error_message = conditionMessage(fit_i)
      )
      next
    }

    out_list[[i]] <- .extract_sensitivity_rows(
      fit = fit_i,
      delta = delta_i,
      method = method,
      baseline_ref = baseline_rows,
      error_message = NA_character_
    )
  }

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  class(out) <- c("spliv_sensitivity_path", "data.frame")

  tipping <- .compute_tipping_points(out)
  attr(out, "call") <- match.call(expand.dots = FALSE)
  attr(out, "delta_grid") <- delta_grid
  attr(out, "method") <- method
  attr(out, "pattern") <- violation_pattern
  attr(out, "tipping_point") <- tipping$values
  attr(out, "tipping_point_message") <- tipping$messages
  out
}

#' Extract Tipping-Point Delta from a Sensitivity Path
#'
#' Returns the tipping-point attribute stored on a
#' `spliv_sensitivity_path()` object. For each term, the tipping point is the
#' smallest `delta` on the supplied grid whose interval includes zero.
#'
#' @param x A `spliv_sensitivity_path` object.
#'
#' @return Named numeric vector of tipping-point values.
#' @examples
#' set.seed(6)
#' d <- data.frame(y = rnorm(80), x = rnorm(80), z = rnorm(80))
#' p <- spliv_sensitivity_path(y ~ x | z, d, method = "uci", delta_grid = c(0, 0.1))
#' spliv_tipping_point(p)
#' @export
spliv_tipping_point <- function(x) {
  if (!inherits(x, "spliv_sensitivity_path")) {
    stop("`x` must inherit from class `spliv_sensitivity_path`.")
  }
  out <- attr(x, "tipping_point")
  attr(out, "message") <- attr(x, "tipping_point_message")
  out
}

.plot_spliv_sensitivity_path <- function(x,
                                         term = NULL,
                                         ylab = "Effect (beta)",
                                         main = NULL,
                                         ...) {
  if (!inherits(x, "spliv_sensitivity_path")) {
    stop("`x` must inherit from class `spliv_sensitivity_path`.")
  }

  terms <- unique(x$term)
  if (is.null(term)) {
    term <- terms[[1]]
    if (length(terms) > 1) {
      warning("Multiple terms are present; plotting the first term `", term, "`.")
    }
  }
  if (!term %in% terms) {
    stop("Requested `term` was not found in the sensitivity path object.")
  }

  dat <- x[x$term == term, , drop = FALSE]
  dat <- dat[order(dat$delta), , drop = FALSE]
  if (is.null(main)) {
    main <- paste0("Sensitivity path: ", term)
  }

  ylim <- range(c(dat$conf_low, dat$conf_high, dat$estimate), na.rm = TRUE)
  graphics::plot(
    dat$delta,
    dat$conf_high,
    type = "n",
    xlab = expression(delta),
    ylab = ylab,
    main = main,
    ylim = ylim,
    ...
  )

  graphics::polygon(
    c(dat$delta, rev(dat$delta)),
    c(dat$conf_low, rev(dat$conf_high)),
    border = NA,
    col = grDevices::adjustcolor("grey60", alpha.f = 0.35)
  )
  graphics::lines(dat$delta, dat$conf_low, lty = 2)
  graphics::lines(dat$delta, dat$conf_high, lty = 2)
  if (any(is.finite(dat$estimate))) {
    graphics::lines(dat$delta, dat$estimate, lwd = 2)
  }
  graphics::abline(h = 0, lty = 3)
  invisible(x)
}

#' Plot Patterned Sensitivity Output or a Fitted Object
#'
#' @param df_plot Data frame returned by a sensitivity helper or `plausexog_fit` object.
#' @param ylab Y-axis label.
#' @param main Plot title.
#' @param ... Unused.
#'
#' @return Invisibly returns the plotted input.
#' @examples
#' set.seed(7)
#' d <- data.frame(y = rnorm(80), x = rnorm(80), z = rnorm(80))
#' p <- spliv_sensitivity_path(y ~ x | z, d, method = "uci", delta_grid = c(0, 0.1))
#' plot_sp_sensitivity(p, term = "x")
#' @export
plot_sp_sensitivity <- function(df_plot,
                                    ylab = "Effect (beta)",
                                    main = "Plausibly Exogenous IV Sensitivity",
                                    ...) {
  if (inherits(df_plot, "spliv_sensitivity_path")) {
    .plot_spliv_sensitivity_path(
      x = df_plot,
      ylab = ylab,
      main = main,
      ...
    )
    return(invisible(df_plot))
  }

  if (inherits(df_plot, "plausexog_fit")) {
    est <- df_plot$estimates
    x <- seq_len(nrow(est))
    graphics::plot(
      x,
      est$estimate %||% rep(NA_real_, nrow(est)),
      xaxt = "n",
      xlab = "Term",
      ylab = ylab,
      main = main,
      pch = 19
    )
    graphics::axis(1, at = x, labels = est$term, las = 2, cex.axis = 0.7)
    if (all(c("conf.low", "conf.high") %in% names(est))) {
      for (i in x) {
        graphics::lines(c(i, i), c(est$conf.low[i], est$conf.high[i]))
      }
    }
    return(invisible(df_plot))
  }

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- ggplot2::ggplot(df_plot, ggplot2::aes(x = delta, group = method)) +
      ggplot2::geom_line(ggplot2::aes(y = conf.low, linetype = method)) +
      ggplot2::geom_line(ggplot2::aes(y = conf.high, linetype = method)) +
      ggplot2::labs(x = expression(delta), y = ylab, title = main, linetype = "Method") +
      ggplot2::theme_minimal()
    if ("estimate" %in% names(df_plot) && any(!is.na(df_plot$estimate))) {
      gg <- gg + ggplot2::geom_line(ggplot2::aes(y = estimate, linetype = method), alpha = 0.6)
    }
    print(gg)
    return(invisible(df_plot))
  }

  methods <- unique(df_plot$method)
  first <- methods[1]
  d1 <- df_plot[df_plot$method == first, , drop = FALSE]
  graphics::plot(
    d1$delta,
    d1$conf.high,
    type = "n",
    xlab = "delta",
    ylab = ylab,
    main = main,
    ylim = range(df_plot$conf.low, df_plot$conf.high, na.rm = TRUE)
  )
  graphics::lines(d1$delta, d1$conf.low)
  graphics::lines(d1$delta, d1$conf.high)
  if ("estimate" %in% names(df_plot) && any(!is.na(d1$estimate))) {
    graphics::lines(d1$delta, d1$estimate)
  }
  if (length(methods) > 1) {
    for (m in methods[-1]) {
      dm <- df_plot[df_plot$method == m, , drop = FALSE]
      graphics::lines(dm$delta, dm$conf.low, lty = 2)
      graphics::lines(dm$delta, dm$conf.high, lty = 2)
      if ("estimate" %in% names(df_plot) && any(!is.na(dm$estimate))) {
        graphics::lines(dm$delta, dm$estimate, lty = 3)
      }
    }
    graphics::legend("topright", legend = methods, lty = 1:length(methods), bty = "n")
  }

  invisible(df_plot)
}

#' @export
plot.plausexog_fit <- function(x, ...) {
  plot_sp_sensitivity(x, ...)
}

#' @export
plot.spliv_sensitivity_path <- function(x,
                                        ylab = "Effect (beta)",
                                        main = "SPLIV Sensitivity Path",
                                        ...) {
  .plot_spliv_sensitivity_path(
    x = x,
    ylab = ylab,
    main = main,
    ...
  )
}

#' @rdname plot_sp_sensitivity
#' @export
plot_conley_sensitivity <- function(df_plot,
                                    ylab = "Effect (beta)",
                                    main = "Plausibly Exogenous IV Sensitivity",
                                    ...) {
  plot_sp_sensitivity(
    df_plot = df_plot,
    ylab = ylab,
    main = main,
    ...
  )
}
