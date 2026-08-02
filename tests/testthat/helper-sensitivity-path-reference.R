reference_sensitivity_path <- function(formula,
                                       data,
                                       method = c("uci", "ltz"),
                                       delta_grid = seq(0, 0.20, by = 0.01),
                                       violation_pattern = NULL,
                                       stop_on_error = TRUE,
                                       ...) {
  method <- match.arg(tolower(method), c("uci", "ltz", "bpe"))
  if (identical(method, "bpe")) {
    stop("Reference sensitivity paths support LTZ and UCI only.")
  }
  delta_grid <- spliv:::.validate_sensitivity_delta_grid(delta_grid)
  dots <- list(...)
  spliv:::.check_sensitivity_path_dots(dots)

  fit_at_delta <- function(delta) {
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

  baseline_fit <- tryCatch(fit_at_delta(0), error = function(e) e)
  if (inherits(baseline_fit, "error")) {
    stop(
      "Failed to fit the baseline sensitivity model at delta = 0: ",
      conditionMessage(baseline_fit),
      call. = FALSE
    )
  }
  baseline_rows <- spliv:::.baseline_reference_rows(baseline_fit, method = method)

  out_list <- vector("list", length(delta_grid))
  zero_tol <- sqrt(.Machine$double.eps)
  for (i in seq_along(delta_grid)) {
    delta_i <- delta_grid[[i]]
    fit_i <- if (abs(delta_i) < zero_tol) {
      baseline_fit
    } else {
      tryCatch(fit_at_delta(delta_i), error = function(e) e)
    }

    if (inherits(fit_i, "error")) {
      if (isTRUE(stop_on_error)) {
        stop(
          "Sensitivity path fit failed at delta = ", format(delta_i),
          ": ", conditionMessage(fit_i),
          call. = FALSE
        )
      }
      out_list[[i]] <- spliv:::.error_sensitivity_rows(
        baseline_ref = baseline_rows,
        delta = delta_i,
        method = method,
        error_message = conditionMessage(fit_i)
      )
    } else {
      out_list[[i]] <- spliv:::.extract_sensitivity_rows(
        fit = fit_i,
        delta = delta_i,
        method = method,
        baseline_ref = baseline_rows,
        error_message = NA_character_
      )
    }
  }

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  class(out) <- c("spliv_sensitivity_path", "data.frame")
  tipping <- spliv:::.compute_tipping_points(out)
  attr(out, "call") <- match.call(expand.dots = FALSE)
  attr(out, "delta_grid") <- delta_grid
  attr(out, "method") <- method
  attr(out, "pattern") <- violation_pattern
  attr(out, "tipping_point") <- tipping$values
  attr(out, "tipping_point_message") <- tipping$messages
  out
}

make_sensitivity_path_optimization_data <- function(n_unit = 18L, n_time = 8L,
                                                    seed = 90210L) {
  set.seed(seed)
  n <- n_unit * n_time
  unit <- factor(rep(seq_len(n_unit), each = n_time))
  time <- factor(rep(seq_len(n_time), times = n_unit))
  unit_effect <- rnorm(n_unit)[unit]
  time_effect <- rnorm(n_time)[time]
  z <- rnorm(n)
  w <- rnorm(n)
  exposure <- stats::pnorm(0.5 * w + rnorm(n, sd = 0.3))
  x <- 0.85 * z + 0.35 * w + unit_effect + time_effect + rnorm(n)
  y <- 0.7 * x + 0.25 * w + 0.12 * exposure * z +
    unit_effect + time_effect + rnorm(n)
  data.frame(y, x, z, w, exposure, unit, time)
}

expect_sensitivity_paths_equivalent <- function(optimized, reference,
                                                tolerance = 1e-10) {
  testthat::expect_identical(class(optimized), class(reference))
  testthat::expect_identical(names(optimized), names(reference))
  testthat::expect_identical(nrow(optimized), nrow(reference))

  numeric_columns <- names(optimized)[vapply(optimized, is.numeric, logical(1))]
  for (column in numeric_columns) {
    testthat::expect_equal(
      optimized[[column]], reference[[column]],
      tolerance = tolerance,
      info = paste("numeric path column", column)
    )
  }
  other_columns <- setdiff(names(optimized), numeric_columns)
  for (column in other_columns) {
    testthat::expect_identical(
      optimized[[column]], reference[[column]],
      info = paste("non-numeric path column", column)
    )
  }

  for (attribute in c(
    "delta_grid", "method", "pattern", "tipping_point",
    "tipping_point_message"
  )) {
    testthat::expect_equal(
      attr(optimized, attribute), attr(reference, attribute),
      tolerance = tolerance,
      info = paste("path attribute", attribute)
    )
  }
  testthat::expect_true(is.call(attr(optimized, "call")))
  invisible(TRUE)
}
