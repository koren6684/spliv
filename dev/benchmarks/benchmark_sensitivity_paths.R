#!/usr/bin/env Rscript

# Fair runtime and allocation benchmark for UCI sensitivity paths. Developer
# material only: dev/ is excluded from the source package.

root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  stop("Run this script from the spliv package root.")
}
for (pkg in c("devtools", "fixest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Install package '", pkg, "' to run the benchmark.")
  }
}
devtools::load_all(root, quiet = TRUE)

seed <- 2400L
fast <- identical(toupper(Sys.getenv("SPLIV_BENCHMARK_FAST", unset = "FALSE")), "TRUE")
iterations <- as.integer(Sys.getenv("SPLIV_BENCHMARK_ITERATIONS", unset = "5"))
theta_steps <- as.integer(Sys.getenv("SPLIV_BENCHMARK_THETA_STEPS", unset = "5"))
if (!is.finite(iterations) || iterations < 5L) {
  stop("SPLIV_BENCHMARK_ITERATIONS must be at least 5.")
}
if (!is.finite(theta_steps) || theta_steps < 2L) {
  stop("SPLIV_BENCHMARK_THETA_STEPS must be at least 2.")
}

case_sizes <- if (fast) {
  data.frame(n = 1000L, groups_1 = 10L, groups_2 = 5L)
} else {
  data.frame(
    n = c(1000L, 1000L, 10000L, 10000L, 50000L, 50000L),
    groups_1 = c(10L, 50L, 50L, 200L, 100L, 500L),
    groups_2 = c(5L, 10L, 10L, 40L, 20L, 100L)
  )
}
grid_lengths <- if (fast) 5L else c(5L, 21L, 41L)
scenarios <- c("uniform", "patterned")
fe_modes <- c("none", "two_way")

make_data <- function(n, groups_1, groups_2, case_seed) {
  set.seed(case_seed)
  group_1 <- factor(sample.int(groups_1, n, replace = TRUE))
  group_2 <- factor(sample.int(groups_2, n, replace = TRUE))
  effect_1 <- rnorm(groups_1)[as.integer(group_1)]
  effect_2 <- rnorm(groups_2)[as.integer(group_2)]
  z <- rnorm(n)
  w <- rnorm(n)
  exposure <- pnorm(w)
  x <- 0.8 * z + 0.3 * w + effect_1 + effect_2 + rnorm(n)
  y <- 1.2 * x + 0.25 * w + 0.15 * exposure * z + effect_1 + effect_2 + rnorm(n)
  data.frame(y, x, z, w, exposure, group_1, group_2)
}

make_pattern <- function(scenario) {
  if (identical(scenario, "uniform")) {
    return(NULL)
  }
  spliv_pattern(
    name = "Benchmark exposure",
    pattern = ~ exposure,
    rationale = "Synthetic benchmark pattern.",
    variables_used = "exposure",
    pattern_type = "theory_defined",
    normalize = "max_abs"
  )
}

reference_setup <- function(d, pattern, fe) {
  baseline <- spliv(
    y ~ x + w | z + w, d, fe = fe, method = "uci", delta = 0,
    vcov = "hc1", scale_instrument = "residual_sd"
  )
  mats <- baseline$internals
  pattern_value <- if (is.null(pattern)) {
    rep(1, nrow(d))
  } else {
    as.numeric(spliv_eval_pattern(pattern, d))
  }
  direct_effect <- pattern_value *
    as.numeric(mats$Z[, "z"]) /
    as.numeric(baseline$residualized_instrument_sd[["z"]])
  processed <- data.frame(
    y = mats$y,
    x = mats$X[, "x"],
    w = mats$X[, "w"],
    z = mats$Z[, "z"],
    direct_effect = direct_effect
  )
  list(
    data = processed,
    formula = if (is.null(fe)) {
      adjusted_y ~ w | 0 | x ~ z
    } else {
      adjusted_y ~ w - 1 | 0 | x ~ z - 1
    }
  )
}

run_spliv <- function(d, pattern, fe, delta_grid) {
  spliv_sensitivity_path(
    y ~ x + w | z + w,
    d,
    fe = fe,
    method = "uci",
    delta_grid = delta_grid,
    violation_pattern = pattern,
    vcov = "hc1",
    scale_instrument = "residual_sd",
    grid = list(steps = theta_steps)
  )
}

run_reference <- function(setup, delta_grid) {
  zcrit <- stats::qnorm(0.975)
  ssc <- fixest::ssc(K.adj = TRUE, K.fixef = "none")
  rows <- lapply(delta_grid, function(delta_i) {
    theta_grid <- seq(-delta_i, delta_i, length.out = theta_steps)
    intervals <- vapply(theta_grid, function(theta_i) {
      benchmark_data <- setup$data
      benchmark_data$adjusted_y <- benchmark_data$y -
        theta_i * benchmark_data$direct_effect
      fit <- fixest::feols(
        setup$formula,
        data = benchmark_data,
        vcov = "hetero",
        ssc = ssc,
        notes = FALSE
      )
      estimate <- unname(stats::coef(fit)[["fit_x"]])
      se <- unname(sqrt(diag(stats::vcov(fit)))[["fit_x"]])
      c(conf_low = estimate - zcrit * se, conf_high = estimate + zcrit * se)
    }, numeric(2))
    data.frame(
      delta = delta_i,
      conf_low = min(intervals["conf_low", ]),
      conf_high = max(intervals["conf_high", ])
    )
  })
  do.call(rbind, rows)
}

allocated_bytes <- function(fun) {
  profile <- tempfile("spliv-rprofmem-", fileext = ".out")
  on.exit(unlink(profile), add = TRUE)
  utils::Rprofmem(profile)
  value <- tryCatch(fun(), finally = utils::Rprofmem(NULL))
  lines <- readLines(profile, warn = FALSE)
  bytes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
  list(value = value, bytes = sum(bytes, na.rm = TRUE))
}

time_repeated <- function(fun, n_iterations) {
  elapsed <- numeric(n_iterations)
  error_text <- character(n_iterations)
  for (i in seq_len(n_iterations)) {
    gc()
    started <- proc.time()[["elapsed"]]
    tryCatch(
      fun(),
      error = function(e) error_text[[i]] <<- conditionMessage(e)
    )
    elapsed[[i]] <- proc.time()[["elapsed"]] - started
  }
  list(elapsed = elapsed, errors = error_text[nzchar(error_text)])
}

summarize_method <- function(method, timing, allocation, metadata, correctness) {
  elapsed <- timing$elapsed
  data.frame(
    metadata,
    method = method,
    iterations = length(elapsed),
    median_seconds = stats::median(elapsed),
    iqr_seconds = stats::IQR(elapsed),
    iterations_per_second = if (stats::median(elapsed) > 0) 1 / stats::median(elapsed) else Inf,
    allocated_bytes = allocation$bytes,
    error_count = length(timing$errors),
    errors = paste(unique(timing$errors), collapse = " | "),
    max_abs_conf_low_difference = correctness[["conf_low"]],
    max_abs_conf_high_difference = correctness[["conf_high"]],
    stringsAsFactors = FALSE
  )
}

message(if (fast) {
  "Running the representative smoke benchmark."
} else {
  "Running the full benchmark matrix; this can take hours on large cases."
})

rows <- list()
counter <- 1L
for (size_i in seq_len(nrow(case_sizes))) {
  for (grid_length in grid_lengths) {
    for (scenario in scenarios) {
      for (fe_mode in fe_modes) {
        n_i <- case_sizes$n[[size_i]]
        g1_i <- case_sizes$groups_1[[size_i]]
        g2_i <- case_sizes$groups_2[[size_i]]
        case_seed <- seed + n_i + g1_i + g2_i + grid_length +
          match(scenario, scenarios) + 10L * match(fe_mode, fe_modes)
        d <- make_data(n_i, g1_i, g2_i, case_seed)
        pattern <- make_pattern(scenario)
        fe <- if (identical(fe_mode, "two_way")) ~ group_1 + group_2 else NULL
        delta_grid <- seq(0, 0.20, length.out = grid_length)
        setup <- reference_setup(d, pattern, fe)
        spliv_fun <- function() run_spliv(d, pattern, fe, delta_grid)
        reference_fun <- function() run_reference(setup, delta_grid)

        message(sprintf(
          "n=%d groups=(%d,%d) delta_grid=%d theta_grid=%d scenario=%s FE=%s",
          n_i, g1_i, g2_i, grid_length, theta_steps, scenario, fe_mode
        ))

        # Warm-up runs are excluded from the timed iterations.
        spliv_warm <- spliv_fun()
        reference_warm <- reference_fun()
        spliv_x <- spliv_warm[spliv_warm$term == "x", c("delta", "conf_low", "conf_high")]
        rownames(spliv_x) <- NULL
        differences <- c(
          conf_low = max(abs(spliv_x$conf_low - reference_warm$conf_low)),
          conf_high = max(abs(spliv_x$conf_high - reference_warm$conf_high))
        )
        if (any(!is.finite(differences)) || any(differences > 1e-8)) {
          stop(
            "The computations are not numerically equivalent for this case. Differences: ",
            paste(names(differences), format(differences), collapse = ", ")
          )
        }

        spliv_timing <- time_repeated(spliv_fun, iterations)
        reference_timing <- time_repeated(reference_fun, iterations)
        spliv_allocation <- allocated_bytes(spliv_fun)
        reference_allocation <- allocated_bytes(reference_fun)
        metadata <- list(
          n = n_i,
          groups_1 = g1_i,
          groups_2 = g2_i,
          fixed_effects = fe_mode,
          delta_grid_length = grid_length,
          theta_grid_length = theta_steps,
          scenario = scenario,
          covariance = "HC1/heteroskedastic",
          treatment_count = 1L,
          instrument_count = 1L,
          control_count = 1L,
          seed = case_seed,
          package_version = as.character(utils::packageVersion("spliv")),
          R_version = R.version.string,
          operating_system = paste(Sys.info()[c("sysname", "release", "machine")], collapse = " ")
        )
        rows[[counter]] <- summarize_method(
          "spliv_sensitivity_path", spliv_timing, spliv_allocation,
          metadata, differences
        )
        counter <- counter + 1L
        rows[[counter]] <- summarize_method(
          "repeated_adjusted_outcome_fixest", reference_timing, reference_allocation,
          metadata, differences
        )
        counter <- counter + 1L
      }
    }
  }
}

out_dir <- file.path(root, "dev", "benchmarks", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, if (fast) "benchmark_results_smoke.csv" else "benchmark_results_full.csv")
utils::write.csv(do.call(rbind, rows), out_file, row.names = FALSE)
message("Wrote ", out_file)
