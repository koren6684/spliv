test_that("optimized paths match repeated public fits across FE, covariance, method, and pattern", {
  d <- make_sensitivity_path_optimization_data()
  formula <- y ~ x + w | z + w
  patterned <- spliv_pattern(
    "Exposure", ~ exposure,
    rationale = "The direct-effect channel follows pre-specified exposure.",
    variables_used = "exposure",
    pattern_type = "theory_defined"
  )

  fe_options <- list(none = NULL, one_way = ~ unit, two_way = ~ unit + time)
  covariance_options <- c("iid", "hc1", "cluster")
  methods <- c("uci", "ltz")
  patterns <- list(uniform = NULL, patterned = patterned)

  for (fe_name in names(fe_options)) {
    for (covariance in covariance_options) {
      for (method in methods) {
        for (pattern_name in names(patterns)) {
          args <- list(
            formula = formula,
            data = d,
            method = method,
            delta_grid = c(0, 0.04, 0.12),
            violation_pattern = patterns[[pattern_name]],
            fe = fe_options[[fe_name]],
            vcov = covariance,
            cluster = if (identical(covariance, "cluster")) ~ unit else NULL,
            grid = list(steps = 5L)
          )
          optimized <- do.call(spliv_sensitivity_path, args)
          reference <- do.call(reference_sensitivity_path, args)
          expect_sensitivity_paths_equivalent(
            optimized, reference, tolerance = 1e-10
          )
        }
      }
    }
  }
})

test_that("complete-case filtering and cluster alignment match repeated fits", {
  base <- make_sensitivity_path_optimization_data(seed = 90211L)
  missing_cases <- list(
    outcome = transform(base, y = replace(y, c(3, 17), NA_real_)),
    treatment = transform(base, x = replace(x, c(5, 19), NA_real_)),
    instrument = transform(base, z = replace(z, c(7, 23), NA_real_))
  )

  for (case_name in names(missing_cases)) {
    d <- missing_cases[[case_name]]
    args <- list(
      formula = y ~ x + w | z + w,
      data = d,
      method = "uci",
      delta_grid = c(0, 0.07),
      vcov = "hc1",
      grid = list(steps = 5L)
    )
    expect_sensitivity_paths_equivalent(
      do.call(spliv_sensitivity_path, args),
      do.call(reference_sensitivity_path, args),
      tolerance = 1e-10
    )
  }

  clustered <- base
  clustered$y[c(2, 31)] <- NA_real_
  clustered$unit[c(11, 57)] <- NA
  cluster_args <- list(
    formula = y ~ x + w | z + w,
    data = clustered,
    method = "ltz",
    delta_grid = c(0, 0.05, 0.1),
    fe = ~ time,
    vcov = "cluster",
    cluster = ~ unit
  )
  optimized <- do.call(spliv_sensitivity_path, cluster_args)
  reference <- do.call(reference_sensitivity_path, cluster_args)
  expect_sensitivity_paths_equivalent(optimized, reference, tolerance = 1e-10)
  expect_true(all(optimized$nobs == nrow(clustered) - 4L))
})

test_that("missing pattern values retain the baseline error behavior", {
  d <- make_sensitivity_path_optimization_data(seed = 90212L)
  d$exposure[9] <- NA_real_
  pattern <- spliv_pattern(
    "Exposure", ~ exposure,
    rationale = "The direct-effect channel follows pre-specified exposure."
  )
  args <- list(
    formula = y ~ x + w | z + w,
    data = d,
    method = "uci",
    delta_grid = c(0, 0.1),
    violation_pattern = pattern
  )

  optimized_error <- tryCatch(do.call(spliv_sensitivity_path, args), error = identity)
  reference_error <- tryCatch(do.call(reference_sensitivity_path, args), error = identity)
  expect_s3_class(optimized_error, "error")
  expect_s3_class(reference_error, "error")
  expect_identical(conditionMessage(optimized_error), conditionMessage(reference_error))
  expect_match(conditionMessage(optimized_error), "contains `NA` values")
})

test_that("formula, function, vector, and column patterns are unchanged", {
  d <- make_sensitivity_path_optimization_data(seed = 90213L)
  patterns <- list(
    formula = spliv_pattern("Formula", ~ exposure, rationale = "Pre-specified."),
    function_pattern = spliv_pattern("Function", function(data) data$exposure,
      rationale = "Pre-specified.", variables_used = "exposure"),
    vector = spliv_pattern("Vector", d$exposure, rationale = "Pre-specified."),
    column = spliv_pattern("Column", "exposure", rationale = "Pre-specified.")
  )

  for (pattern in patterns) {
    args <- list(
      formula = y ~ x + w | z + w,
      data = d,
      method = "uci",
      delta_grid = c(0, 0.08),
      violation_pattern = pattern,
      grid = list(steps = 5L)
    )
    expect_sensitivity_paths_equivalent(
      do.call(spliv_sensitivity_path, args),
      do.call(reference_sensitivity_path, args),
      tolerance = 1e-10
    )
  }
})

test_that("residual-SD paths are unchanged and invariant to instrument units", {
  d <- make_sensitivity_path_optimization_data(seed = 90214L)
  results <- lapply(c(1, 100, 0.01), function(multiplier) {
    scaled <- d
    scaled$z_scaled <- multiplier * scaled$z
    args <- list(
      formula = y ~ x + w | z_scaled + w,
      data = scaled,
      method = "uci",
      delta_grid = c(0, 0.03, 0.09),
      scale_instrument = "residual_sd",
      vcov = "hc1",
      grid = list(steps = 7L)
    )
    optimized <- do.call(spliv_sensitivity_path, args)
    reference <- do.call(reference_sensitivity_path, args)
    expect_sensitivity_paths_equivalent(optimized, reference, tolerance = 1e-10)
    optimized
  })

  for (candidate in results[-1L]) {
    expect_equal(candidate$estimate, results[[1L]]$estimate, tolerance = 1e-10)
    expect_equal(candidate$conf_low, results[[1L]]$conf_low, tolerance = 1e-10)
    expect_equal(candidate$conf_high, results[[1L]]$conf_high, tolerance = 1e-10)
    expect_identical(candidate$contains_zero, results[[1L]]$contains_zero)
  }

  raw_args <- list(
    formula = y ~ x + w | z + w,
    data = d,
    method = "ltz",
    delta_grid = c(0, 0.06),
    scale_instrument = "none",
    vcov = "iid"
  )
  expect_sensitivity_paths_equivalent(
    do.call(spliv_sensitivity_path, raw_args),
    do.call(reference_sensitivity_path, raw_args),
    tolerance = 1e-10
  )
})

test_that("short and long grids preserve values, metadata, and tipping points", {
  d <- make_sensitivity_path_optimization_data(seed = 90215L)
  for (delta_grid in list(c(0.05), seq(0, 0.24, length.out = 13L))) {
    args <- list(
      formula = y ~ x + w | z + w,
      data = d,
      method = "uci",
      delta_grid = delta_grid,
      vcov = "hc1",
      grid = list(steps = 7L)
    )
    optimized <- do.call(spliv_sensitivity_path, args)
    reference <- do.call(reference_sensitivity_path, args)
    expect_sensitivity_paths_equivalent(optimized, reference, tolerance = 1e-10)
    expect_equal(spliv_tipping_point(optimized), spliv_tipping_point(reference),
      tolerance = 1e-10)
    expect_identical(attr(optimized, "delta_grid"), as.numeric(delta_grid))
    expect_identical(attr(optimized, "method"), "uci")
    expect_null(attr(optimized, "pattern"))
  }
})

test_that("formula, FE, cluster, pattern, and IV preparation occur once per path", {
  d <- make_sensitivity_path_optimization_data(seed = 90216L)
  pattern <- spliv_pattern(
    "Exposure", ~ exposure,
    rationale = "The direct-effect channel follows pre-specified exposure."
  )
  counts <- new.env(parent = emptyenv())
  hook <- function(event) {
    current <- if (exists(event, envir = counts, inherits = FALSE)) {
      get(event, envir = counts, inherits = FALSE)
    } else {
      0L
    }
    assign(event, current + 1L, envir = counts)
  }
  old_options <- options(spliv.sensitivity_path.preparation_hook = hook)
  on.exit(options(old_options), add = TRUE)

  spliv_sensitivity_path(
    y ~ x + w | z + w, d,
    method = "uci",
    delta_grid = seq(0, 0.2, length.out = 9L),
    violation_pattern = pattern,
    fe = ~ unit + time,
    vcov = "cluster",
    cluster = ~ unit,
    grid = list(steps = 5L)
  )

  expected_events <- c(
    "formula_preparation", "complete_case_filtering",
    "fixed_effect_residualization", "cluster_alignment",
    "instrument_scaling_preparation", "pattern_preparation",
    "iv_design_preparation"
  )
  expect_identical(sort(ls(counts)), sort(expected_events))
  for (event in expected_events) {
    expect_identical(get(event, envir = counts), 1L, info = event)
  }
})
