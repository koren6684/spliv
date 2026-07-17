.get_cluster_id <- function(cluster, data_cc) {
  if (is.null(cluster)) {
    return(NULL)
  }
  if (inherits(cluster, "formula")) {
    mf <- stats::model.frame(cluster, data = data_cc, na.action = stats::na.pass)
    if (ncol(mf) != 1) {
      stop("cluster formula must resolve to exactly one variable, e.g. ~ gid")
    }
    return(mf[[1]])
  }
  if (length(cluster) != nrow(data_cc)) {
    stop("cluster vector length must match number of complete-case rows used in estimation.")
  }
  cluster
}

.build_W_from_XZ <- function(X, Z) {
  x_names <- colnames(X)
  z_names <- colnames(Z)
  common <- intersect(x_names, z_names)
  common <- setdiff(common, "(Intercept)")
  if (length(common) == 0) {
    return(matrix(0, nrow = nrow(X), ncol = 0))
  }
  X[, common, drop = FALSE]
}

.default_uci_inst <- function(X, Z) {
  z_names <- colnames(Z)
  x_names <- colnames(X)
  cand <- setdiff(z_names, c(x_names, "(Intercept)"))
  if (length(cand) == 0) {
    cand <- setdiff(z_names, "(Intercept)")
  }
  cand[1]
}

.first_stage_target <- function(X, Z) {
  x_names <- setdiff(colnames(X), "(Intercept)")
  z_names <- setdiff(colnames(Z), "(Intercept)")
  cand <- setdiff(x_names, z_names)
  if (length(cand) > 0) {
    return(cand[1])
  }
  if (length(x_names) == 0) {
    stop("Could not determine a first-stage target from X.")
  }
  x_names[1]
}

.make_na_ltz_table <- function(coef_names) {
  data.frame(
    term = coef_names,
    estimate = NA_real_,
    std.error = NA_real_,
    conf.low = NA_real_,
    conf.high = NA_real_,
    row.names = NULL
  )
}

.make_not_applicable_fit <- function(call,
                                     formula,
                                     fe,
                                     fe_engine,
                                     vcov,
                                     cluster,
                                     method_out,
                                     method_core,
                                     coef_names,
                                     prior,
                                     grid,
                                     bpe_diag,
                                     msg,
                                     y,
                                     X,
                                     Z,
                                     W,
                                     cluster_id,
                                     scale_instrument,
                                     residualized_instrument_sd) {
  k <- length(coef_names)
  fit <- list(
    call = call,
    formula = formula,
    fe = fe,
    fe_engine = if (is.null(fe)) NULL else fe_engine,
    vcov = vcov,
    cluster = cluster,
    method = method_out,
    method_core = method_core,
    estimates = .make_na_ltz_table(coef_names),
    beta_hat = rep(NA_real_, k),
    beta_iv = rep(NA_real_, k),
    vcov_beta = matrix(NA_real_, k, k),
    prior = prior,
    mu_used = prior$mu %||% NULL,
    Omega_used = prior$Omega %||% NULL,
    A = NULL,
    A_mu = NULL,
    beta_shift = NULL,
    grid = grid,
    bpe = TRUE,
    diag = list(bpe = bpe_diag),
    bpe_diagnostics = bpe_diag,
    scale_instrument = scale_instrument,
    residualized_instrument_sd = residualized_instrument_sd,
    errors = msg,
    warnings = msg,
    internals = list(y = y, X = X, Z = Z, W = W, cluster_id = cluster_id)
  )
  class(fit) <- "plausexog_fit"
  fit
}

.align_bpe_prior_to_Z <- function(validation, Z) {
  z_current <- colnames(Z)
  if (is.null(z_current)) {
    z_current <- paste0("z", seq_len(ncol(Z)))
    colnames(Z) <- z_current
  }

  mu_full <- validation$prior_mu_full
  Omega_full <- validation$prior_Omega_full
  if (!is.null(mu_full) &&
      !is.null(names(mu_full)) &&
      is.matrix(Omega_full) &&
      !is.null(rownames(Omega_full)) &&
      !is.null(colnames(Omega_full)) &&
      all(z_current %in% names(mu_full)) &&
      all(z_current %in% rownames(Omega_full)) &&
      all(z_current %in% colnames(Omega_full))) {
    mu <- as.numeric(mu_full[z_current])
    names(mu) <- z_current
    Omega <- as.matrix(Omega_full[z_current, z_current, drop = FALSE])
    return(list(mu = mu, Omega = Omega))
  }

  mu_sub <- validation$prior_mu_sub
  Omega_sub <- validation$prior_Omega_sub
  if (length(z_current) == 1L &&
      length(mu_sub) == 1L &&
      is.matrix(Omega_sub) &&
      all(dim(Omega_sub) == c(1L, 1L))) {
    mu <- stats::setNames(as.numeric(mu_sub), z_current)
    Omega <- matrix(
      as.numeric(Omega_sub[1, 1]),
      nrow = 1L,
      ncol = 1L,
      dimnames = list(z_current, z_current)
    )
    return(list(mu = mu, Omega = Omega))
  }

  stop(
    "Internal BPE prior alignment failed: the validated BPE direct-effect prior ",
    "could not be matched to the residualized instrument columns used for estimation."
  )
}

.coerce_bpe_design <- function(bpe_design_arg, bpe_spec, parsed, data) {
  bpe_spec <- bpe_spec %||% list()
  if (!is.list(bpe_spec)) {
    stop("`bpe_spec` must be a list.")
  }

  if (!is.null(bpe_spec$subset_rule)) {
    stop(
      paste(
        "Confirmatory BPE no longer accepts `bpe_spec$subset_rule`.",
        "Use `bpe_explore_subsets()` for exploratory work, then create a",
        "theory-justified `bpe_design()` object for confirmatory estimation."
      )
    )
  }
  if (!is.null(bpe_spec$find_subset) || !is.null(bpe_spec$f_threshold)) {
    stop(
      "Deprecated exploratory subset-search controls are not accepted for confirmatory BPE. Use `bpe_design()` instead."
    )
  }

  design_from_spec <- bpe_spec$design %||% NULL
  raw_subset <- bpe_spec$subset %||% NULL

  n_designs <- sum(!vapply(list(bpe_design_arg, design_from_spec, raw_subset), is.null, logical(1)))
  if (n_designs > 1) {
    stop(
      "Confirmatory BPE requires exactly one declared design/subset source: ",
      "`bpe_design`, `bpe_spec$design`, or `bpe_spec$subset`."
    )
  }

  design_obj <- bpe_design_arg %||% design_from_spec
  if (!is.null(design_obj)) {
    if (!inherits(design_obj, "spliv_bpe_design")) {
      stop("`bpe_design` and `bpe_spec$design` must inherit from class `spliv_bpe_design`.")
    }
    return(design_obj)
  }

  if (is.null(raw_subset)) {
    stop(
      paste(
        "Confirmatory BPE requires a researcher-supplied `bpe_design()` object",
        "or an explicit subset supplied through `bpe_spec$subset` together with",
        "confirmatory metadata such as `rationale` and `pre_specified`."
      )
    )
  }

  if (is.logical(raw_subset) && length(raw_subset) == length(parsed$y)) {
    full_subset <- rep(FALSE, nrow(data))
    full_subset[parsed$keep] <- raw_subset
    raw_subset <- full_subset
  }

  rationale <- bpe_spec$rationale %||% NULL
  if (!is.character(rationale) || length(rationale) != 1 || !nzchar(trimws(rationale))) {
    stop("Raw `bpe_spec$subset` requires a non-empty `rationale` for confirmatory BPE.")
  }
  if (is.null(bpe_spec$pre_specified)) {
    stop("Raw `bpe_spec$subset` requires explicit `pre_specified = TRUE` for confirmatory BPE.")
  }
  if (!isTRUE(bpe_spec$pre_specified)) {
    stop("Raw `bpe_spec$subset` requires `pre_specified = TRUE` for confirmatory BPE.")
  }

  bpe_design(
    name = bpe_spec$name %||% "Researcher-supplied BPE subset",
    subset = raw_subset,
    rationale = rationale,
    variables_used = bpe_spec$variables_used %||% NULL,
    subset_type = bpe_spec$subset_type %||% NULL,
    pre_specified = bpe_spec$pre_specified,
    transportability_rationale = bpe_spec$transportability_rationale %||% NULL,
    notes = bpe_spec$notes %||% NULL
  )
}

.instrument_residual_sds <- function(X, Z, W, fe_present) {
  inst_names <- setdiff(colnames(Z), c(colnames(X), "(Intercept)"))
  if (length(inst_names) == 0) {
    inst_names <- setdiff(colnames(Z), "(Intercept)")
  }
  if (length(inst_names) == 0) {
    return(setNames(numeric(0), character(0)))
  }

  out <- vapply(
    inst_names,
    function(z_name) {
      .bpe_residualized_sd(
        Z_sub = Z,
        W_sub = W,
        z_name = z_name,
        fe_present = fe_present
      )
    },
    numeric(1)
  )
  setNames(as.numeric(out), inst_names)
}

.delta_bounds_from_scale <- function(delta, inst, residualized_instrument_sd, scale_instrument) {
  p <- length(inst)
  if (length(delta) == 1) {
    delta <- rep(delta, p)
  }
  if (length(delta) != p) {
    stop("`grid$delta` must be scalar or have one value per varied instrument.")
  }
  delta <- as.numeric(delta)
  if (any(!is.finite(delta)) || any(delta < 0)) {
    stop("`grid$delta` must contain non-negative finite values.")
  }

  if (identical(scale_instrument, "none")) {
    gamma_abs <- delta
  } else {
    sd_inst <- residualized_instrument_sd[inst]
    if (any(!is.finite(sd_inst)) || any(sd_inst <= 0)) {
      stop(
        "Residualized instrument SD is required and must be positive when `scale_instrument = 'residual_sd'`."
      )
    }
    gamma_abs <- delta / sd_inst
  }

  list(gmin = -gamma_abs, gmax = gamma_abs)
}

#' Create a Patterned Exclusion-Violation Object
#'
#' Creates a researcher-specified direct-effect pattern for patterned
#' sensitivity analysis. The pattern determines where direct effects of the
#' instrument are expected to be larger or smaller in the estimation sample.
#'
#' @param name Short pattern name.
#' @param pattern Pattern specification. Accepts a one-sided formula, a
#'   `function(data)`, a numeric vector of length `nrow(data)`, a logical vector
#'   of length `nrow(data)`, or a character string naming a numeric or logical
#'   column in `data`.
#' @param rationale Substantive justification for the proposed pattern.
#' @param variables_used Optional character vector describing the variables used
#'   to define the pattern.
#' @param pattern_type Optional short label such as `"theory_defined"` or
#'   `"design_based"`.
#' @param normalize One of `"max_abs"` (default), `"sd"`, or `"none"`.
#' @param center Logical; if `TRUE`, subtracts the estimation-sample mean before
#'   normalization.
#' @param pre_specified Logical indicator for whether the pattern was chosen
#'   before examining outcome results.
#' @param notes Optional free-form notes.
#'
#' @return An object of class `"spliv_pattern"`.
#' @examples
#' p <- spliv_pattern(
#'   name = "Exposure", pattern = ~ exposure,
#'   rationale = "The alternative channel follows exposure.",
#'   variables_used = "exposure", pattern_type = "theory_defined"
#' )
#' p$name
#' @export
spliv_pattern <- function(name,
                          pattern,
                          rationale,
                          variables_used = NULL,
                          pattern_type = NULL,
                          normalize = c("max_abs", "sd", "none"),
                          center = FALSE,
                          pre_specified = TRUE,
                          notes = NULL) {
  if (!is.character(name) || length(name) != 1 || !nzchar(trimws(name))) {
    stop("`name` must be a non-empty character scalar.")
  }
  if (missing(pattern)) {
    stop("`pattern` must be supplied.")
  }
  if (!is.character(rationale) || length(rationale) != 1 || !nzchar(trimws(rationale))) {
    stop("`rationale` must be a non-empty character scalar.")
  }
  if (!is.logical(center) || length(center) != 1 || is.na(center)) {
    stop("`center` must be TRUE or FALSE.")
  }
  if (!is.logical(pre_specified) || length(pre_specified) != 1 || is.na(pre_specified)) {
    stop("`pre_specified` must be TRUE or FALSE.")
  }

  normalize <- match.arg(normalize)
  pattern_kind <- .spliv_pattern_kind(pattern)
  if (is.null(variables_used)) {
    variables_used <- .spliv_pattern_infer_variables(pattern, pattern_kind = pattern_kind)
  }
  if (!is.null(variables_used)) {
    variables_used <- unique(as.character(variables_used))
    variables_used <- variables_used[nzchar(trimws(variables_used))]
    if (!length(variables_used)) {
      variables_used <- NULL
    }
  }

  out <- list(
    name = trimws(name),
    pattern = pattern,
    rationale = rationale,
    variables_used = variables_used,
    pattern_type = pattern_type,
    normalize = normalize,
    center = center,
    pre_specified = pre_specified,
    notes = notes,
    pattern_kind = pattern_kind,
    created_at = Sys.time(),
    call = match.call()
  )
  class(out) <- "spliv_pattern"
  out
}

.spliv_pattern_kind <- function(pattern) {
  if (inherits(pattern, "formula")) {
    if (length(pattern) != 2) {
      stop("`pattern` formulas must be one-sided, e.g. `~ treatment_channel_exposure`.")
    }
    return("formula")
  }
  if (is.function(pattern)) {
    return("function")
  }
  if (is.numeric(pattern)) {
    return("numeric")
  }
  if (is.logical(pattern)) {
    return("logical")
  }
  if (is.character(pattern) && length(pattern) == 1 && nzchar(trimws(pattern))) {
    return("column")
  }
  stop(
    "`pattern` must be a one-sided formula, a function(data), a numeric vector, ",
    "a logical vector, or a character string naming a numeric/logical column."
  )
}

.spliv_pattern_infer_variables <- function(pattern, pattern_kind = .spliv_pattern_kind(pattern)) {
  if (identical(pattern_kind, "formula")) {
    return(all.vars(pattern))
  }
  if (identical(pattern_kind, "column")) {
    return(as.character(pattern))
  }
  NULL
}

.spliv_pattern_variables <- function(pattern) {
  vars <- pattern$variables_used %||% .spliv_pattern_infer_variables(
    pattern$pattern,
    pattern_kind = pattern$pattern_kind
  )
  if (is.null(vars)) {
    return(character(0))
  }
  unique(as.character(vars))
}

.spliv_eval_formula_pattern <- function(formula, data) {
  eval(formula[[2L]], envir = data, enclos = parent.frame())
}

.spliv_eval_function_pattern <- function(fun, data) {
  fun(data)
}

.spliv_eval_column_pattern <- function(col_name, data) {
  if (!col_name %in% names(data)) {
    stop("Pattern column `", col_name, "` was not found in `data`.")
  }
  data[[col_name]]
}

.spliv_pattern_summary <- function(x) {
  x <- as.numeric(x)
  list(
    n = length(x),
    mean = mean(x),
    sd = if (length(x) <= 1) 0 else stats::sd(x),
    min = min(x),
    max = max(x)
  )
}

.spliv_eval_pattern_sample <- function(pattern, data, sample_idx = NULL) {
  if (!inherits(pattern, "spliv_pattern")) {
    stop("`pattern` must inherit from class `spliv_pattern`.")
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  }
  if (is.null(sample_idx)) {
    sample_idx <- rep(TRUE, nrow(data))
  }
  if (!is.logical(sample_idx) || length(sample_idx) != nrow(data)) {
    stop("`sample_idx` must be a logical vector of length `nrow(data)`.")
  }
  if (!any(sample_idx)) {
    stop("`sample_idx` selects no observations.")
  }

  raw_full <- tryCatch(
    {
      if (identical(pattern$pattern_kind, "formula")) {
        .spliv_eval_formula_pattern(pattern$pattern, data)
      } else if (identical(pattern$pattern_kind, "function")) {
        .spliv_eval_function_pattern(pattern$pattern, data)
      } else if (identical(pattern$pattern_kind, "numeric")) {
        pattern$pattern
      } else if (identical(pattern$pattern_kind, "logical")) {
        pattern$pattern
      } else if (identical(pattern$pattern_kind, "column")) {
        .spliv_eval_column_pattern(pattern$pattern, data)
      } else {
        stop("Unsupported pattern kind: ", pattern$pattern_kind)
      }
    },
    error = function(e) {
      stop("Failed to evaluate violation pattern `", pattern$name, "`: ", conditionMessage(e), call. = FALSE)
    }
  )

  allow_scalar_recycle <- identical(pattern$pattern_kind, "formula") || identical(pattern$pattern_kind, "function")
  if (length(raw_full) == 1 && allow_scalar_recycle) {
    raw_full <- rep(raw_full, nrow(data))
  }
  if (length(raw_full) != nrow(data)) {
    stop(
      "Violation pattern `", pattern$name, "` returned length ", length(raw_full),
      " but `data` has ", nrow(data), " rows."
    )
  }
  if (!(is.numeric(raw_full) || is.logical(raw_full))) {
    stop("Violation pattern `", pattern$name, "` must evaluate to a numeric or logical vector.")
  }

  raw_full <- if (is.logical(raw_full)) as.numeric(raw_full) else as.numeric(raw_full)
  raw_sample <- raw_full[sample_idx]

  if (all(is.na(raw_sample))) {
    stop("Violation pattern `", pattern$name, "` evaluated to all `NA` on the estimation sample.")
  }
  if (any(is.na(raw_sample))) {
    stop("Violation pattern `", pattern$name, "` contains `NA` values on the estimation sample.")
  }
  if (any(is.infinite(raw_sample))) {
    stop("Violation pattern `", pattern$name, "` contains infinite values on the estimation sample.")
  }
  if (all(raw_sample == 0)) {
    stop("Violation pattern `", pattern$name, "` is all zero on the estimation sample.")
  }

  work_sample <- raw_sample
  if (isTRUE(pattern$center)) {
    work_sample <- work_sample - mean(work_sample)
  }
  if (all(work_sample == 0)) {
    stop("Violation pattern `", pattern$name, "` becomes all zero after centering on the estimation sample.")
  }

  warnings <- character(0)
  scale_factor <- 1
  if (identical(pattern$normalize, "max_abs")) {
    scale_factor <- max(abs(work_sample))
    if (!is.finite(scale_factor) || scale_factor <= 0) {
      stop("Violation pattern `", pattern$name, "` cannot be normalized by `max_abs`.")
    }
  } else if (identical(pattern$normalize, "sd")) {
    scale_factor <- stats::sd(work_sample)
    if (!is.finite(scale_factor) || scale_factor <= 0) {
      stop("Violation pattern `", pattern$name, "` cannot be normalized by `sd` because the sample SD is zero.")
    }
  } else if (identical(pattern$normalize, "none")) {
    warnings <- c(
      warnings,
      paste0(
        "Violation pattern `", pattern$name, "` uses `normalize = 'none'`; ",
        "theta/delta depends on the supplied pattern scale."
      )
    )
  }

  norm_sample <- if (identical(pattern$normalize, "none")) work_sample else work_sample / scale_factor

  list(
    values = as.numeric(norm_sample),
    raw_values = as.numeric(raw_sample),
    raw_summary = .spliv_pattern_summary(raw_sample),
    normalized_summary = .spliv_pattern_summary(norm_sample),
    warnings = warnings,
    sample_n = length(norm_sample),
    variables_used = .spliv_pattern_variables(pattern)
  )
}

#' Evaluate a Direct-Effect Pattern
#'
#' Evaluates a `spliv_pattern()` object on a data frame and returns a numeric
#' vector after optional centering and normalization.
#'
#' @param pattern A `spliv_pattern()` object.
#' @param data Data frame used for evaluation.
#'
#' @return Numeric vector of length `nrow(data)`.
#' @examples
#' d <- data.frame(exposure = seq(-1, 1, length.out = 5))
#' p <- spliv_pattern("Exposure", ~ exposure,
#'   rationale = "Illustration of a monotone exposure pattern.")
#' spliv_eval_pattern(p, d)
#' @export
spliv_eval_pattern <- function(pattern, data) {
  out <- .spliv_eval_pattern_sample(pattern = pattern, data = data, sample_idx = rep(TRUE, nrow(data)))
  structure(
    out$values,
    spliv_pattern_raw_summary = out$raw_summary,
    spliv_pattern_summary = out$normalized_summary,
    spliv_pattern_warnings = out$warnings
  )
}

.spliv_excluded_instruments <- function(X, Z) {
  z_names <- setdiff(colnames(Z), "(Intercept)")
  x_names <- setdiff(colnames(X), "(Intercept)")
  excluded <- setdiff(z_names, x_names)
  if (!length(excluded)) {
    excluded <- z_names
  }
  excluded
}

.spliv_pattern_targets <- function(X, Z) {
  x_names <- setdiff(colnames(X), "(Intercept)")
  z_names <- .spliv_excluded_instruments(X, Z)

  target_x <- setdiff(x_names, setdiff(colnames(Z), "(Intercept)"))
  if (length(target_x) != 1) {
    stop(
      "Patterned sensitivity currently supports one endogenous treatment and one instrument. ",
      "Please fit a model with exactly one endogenous regressor."
    )
  }
  if (length(z_names) != 1) {
    stop(
      "Patterned sensitivity currently supports one endogenous treatment and one instrument. ",
      "Please fit a model with exactly one excluded instrument."
    )
  }

  list(x_name = target_x[[1]], z_name = z_names[[1]])
}

.spliv_pattern_scale_vector <- function(z_vec, z_name, residualized_instrument_sd, scale_instrument) {
  z_vec <- as.numeric(z_vec)
  if (identical(scale_instrument, "none")) {
    return(list(values = z_vec, scale_factor = 1))
  }
  sd_inst <- residualized_instrument_sd[[z_name]]
  if (!is.finite(sd_inst) || sd_inst <= 0) {
    stop(
      "Residualized instrument SD is required and must be positive when `scale_instrument = 'residual_sd'`."
    )
  }
  list(values = z_vec / sd_inst, scale_factor = sd_inst)
}

.coerce_pattern_ltz_prior <- function(prior, delta = NULL) {
  if (!is.null(prior) && !is.null(delta)) {
    stop("Supply either `prior` or `delta` for LTZ patterned sensitivity, not both.")
  }

  if (!is.null(prior)) {
    mu <- prior$mu %||% prior$mu_hat %||% NULL
    Omega <- prior$Omega %||% prior$omega %||% prior$omega_hat %||% NULL
    if (is.null(mu) || is.null(Omega)) {
      stop("Patterned LTZ requires scalar `prior$mu` and `prior$Omega` (or `omega`).")
    }
    mu <- as.numeric(mu)
    if (length(mu) != 1 || !is.finite(mu)) {
      stop("Patterned LTZ currently requires a scalar prior mean for the violation pattern.")
    }
    if (is.numeric(Omega) && length(Omega) == 1) {
      Omega <- matrix(as.numeric(Omega), 1, 1)
    }
    if (!is.matrix(Omega) || any(dim(Omega) != c(1, 1)) || !is.finite(Omega[1, 1]) || Omega[1, 1] < 0) {
      stop("Patterned LTZ currently requires a 1 x 1 non-negative prior covariance for the violation pattern.")
    }
    return(list(mu = mu, Omega = Omega))
  }

  if (is.null(delta)) {
    stop("Patterned LTZ requires either `delta` or a scalar `prior`.")
  }
  delta <- as.numeric(delta)
  if (length(delta) != 1 || !is.finite(delta) || delta < 0) {
    stop("`delta` must be a non-negative finite scalar for patterned LTZ.")
  }
  list(mu = 0, Omega = matrix(delta^2, 1, 1))
}

.scalar_ltz_prior_from_delta <- function(delta, X, Z, residualized_instrument_sd, scale_instrument) {
  delta <- as.numeric(delta)
  if (length(delta) != 1 || !is.finite(delta) || delta < 0) {
    stop("`delta` must be a non-negative finite scalar.")
  }
  excluded <- .spliv_excluded_instruments(X, Z)
  if (length(excluded) != 1) {
    stop(
      "Automatic scalar LTZ sensitivity with `delta` currently requires exactly one excluded instrument. ",
      "Otherwise supply an explicit `prior`."
    )
  }

  scale_sd <- if (identical(scale_instrument, "none")) {
    1
  } else {
    sd_inst <- residualized_instrument_sd[[excluded[[1]]]]
    if (!is.finite(sd_inst) || sd_inst <= 0) {
      stop(
        "Residualized instrument SD is required and must be positive when `scale_instrument = 'residual_sd'`."
      )
    }
    sd_inst
  }

  mu <- setNames(rep(0, ncol(Z)), colnames(Z))
  Omega <- matrix(0, ncol(Z), ncol(Z), dimnames = list(colnames(Z), colnames(Z)))
  Omega[excluded[[1]], excluded[[1]]] <- (delta / scale_sd)^2
  list(mu = mu, Omega = Omega, inst = excluded[[1]])
}

.prepare_violation_pattern <- function(violation_pattern,
                                       data,
                                       sample_idx,
                                       X,
                                       Z,
                                       residualized_instrument_sd,
                                       scale_instrument) {
  if (!inherits(violation_pattern, "spliv_pattern")) {
    stop("`violation_pattern` must inherit from class `spliv_pattern`.")
  }

  targets <- .spliv_pattern_targets(X, Z)
  pattern_eval <- .spliv_eval_pattern_sample(
    pattern = violation_pattern,
    data = data,
    sample_idx = sample_idx
  )
  z_scaled <- .spliv_pattern_scale_vector(
    z_vec = Z[, targets$z_name, drop = TRUE],
    z_name = targets$z_name,
    residualized_instrument_sd = residualized_instrument_sd,
    scale_instrument = scale_instrument
  )

  direct_effect <- as.numeric(pattern_eval$values) * as.numeric(z_scaled$values)
  direct_effect_mat <- matrix(
    direct_effect,
    ncol = 1,
    dimnames = list(NULL, paste0("theta_", targets$z_name, "_pattern"))
  )

  info <- list(
    name = violation_pattern$name,
    rationale = violation_pattern$rationale,
    variables_used = pattern_eval$variables_used,
    pattern_type = violation_pattern$pattern_type,
    normalize = violation_pattern$normalize,
    center = violation_pattern$center,
    pre_specified = violation_pattern$pre_specified,
    notes = violation_pattern$notes,
    sample_n = pattern_eval$sample_n,
    raw_pattern_summary = pattern_eval$raw_summary,
    normalized_pattern_summary = pattern_eval$normalized_summary,
    direct_effect_regressor_summary = .spliv_pattern_summary(direct_effect),
    delta_theta_interpretation = if (identical(scale_instrument, "residual_sd")) {
      paste(
        "Theta/Delta is the outcome-unit direct effect for a one residual-standard-deviation",
        "shift in the instrument, scaled by the supplied normalized pattern."
      )
    } else {
      paste(
        "Theta/Delta is the outcome-unit direct effect per unit of the instrument column",
        "used in the estimating equations, scaled by the supplied normalized pattern."
      )
    },
    residual_sd_scaling_used = identical(scale_instrument, "residual_sd"),
    instrument = targets$z_name,
    endogenous_treatment = targets$x_name,
    warnings = pattern_eval$warnings
  )

  list(
    direct_effect = direct_effect_mat,
    direct_effect_names = colnames(direct_effect_mat),
    info = info,
    instrument = targets$z_name,
    endogenous_treatment = targets$x_name,
    warnings = pattern_eval$warnings
  )
}

#' Patterned Sensitivity Analysis for Plausibly Exogenous IV
#'
#' Main estimator for patterned sensitivity analysis of exclusion violations in
#' fixed-effect or residualized IV designs.
#'
#' @param formula IV formula `y ~ X | Z`.
#' @param data Data frame.
#' @param fe Optional one-sided formula of fixed effects.
#' @param fe_engine FE demeaning engine, one of `"fixest"` or `"lfe"`.
#' @param vcov One of `"iid"`, `"hc1"`, or `"cluster"`.
#' @param cluster Cluster ids or one-sided formula, required when `vcov = "cluster"`.
#' @param method One of `"ltz"`, `"uci"`, or `"bpe"`.
#' @param prior Optional prior list with `mu` and `Omega` (or `omega`) for LTZ.
#'   When `violation_pattern` is supplied, patterned LTZ currently requires a
#'   scalar prior over the pattern coefficient.
#' @param delta Optional non-negative scalar sensitivity magnitude. With
#'   `method = "uci"`, `delta` implies theta bounds `[-delta, +delta]` unless
#'   explicit bounds are supplied in `grid`. With `method = "ltz"` and no
#'   explicit `prior`, `delta` induces a zero-mean normal LTZ prior.
#' @param violation_pattern Optional `spliv_pattern()` object describing how the
#'   direct effect of the instrument may vary across observations. If omitted,
#'   LTZ/UCI retain the package's backward-compatible uniform direct-effect
#'   behavior. This argument is currently supported for LTZ and UCI, but not for
#'   confirmatory BPE.
#' @param bpe Logical; when `TRUE`, learn prior moments from a confirmatory
#'   `bpe_design()` and run LTZ.
#' @param bpe_design Optional `bpe_design()` object for confirmatory BPE.
#' @param bpe_spec Optional list. Prefer `bpe_spec = list(design = my_design)`.
#'   For backward compatibility, `bpe_spec$subset` may also be supplied, but it
#'   must represent an explicit researcher-supplied subset and should be paired
#'   with a non-empty `rationale`, explicit `pre_specified = TRUE`, and any
#'   relevant metadata such as `variables_used` or `subset_type`. Supply exactly
#'   one confirmatory subset source: `bpe_design`, `bpe_spec$design`, or
#'   `bpe_spec$subset`. Exploratory `subset_rule` inputs are not accepted for
#'   confirmatory estimation.
#' @param bpe_kappa Positive scalar multiplier applied to the confirmatory BPE
#'   covariance after transport adjustment.
#' @param bpe_omega Deprecated. Confirmatory BPE now uses the full reduced-form
#'   covariance from the subset together with `bpe_transport`.
#' @param bpe_min_n_S Minimum subset size required for BPE eligibility. Default
#'   `2000`.
#' @param bpe_min_clusters_S Minimum number of clusters required in subset `S`
#'   when `vcov = "cluster"`. Default `30`.
#' @param bpe_max_F_S Deprecated. The first-stage F-statistic is reported for
#'   diagnostics only and no longer determines confirmatory BPE eligibility.
#' @param bpe_min_varZ_S Minimum residualized instrument variance required in
#'   subset `S`. Default `1e-6`.
#' @param bpe_equiv_margin Researcher-specified first-stage equivalence margin.
#'   Confirmatory BPE currently supports exactly one instrument.
#' @param bpe_equiv_level Confidence level for the first-stage equivalence check.
#' @param bpe_transport One of `"none"`, `"sampling"`, or `"conservative"`.
#' @param bpe_transport_kappa Non-negative scalar controlling the conservative
#'   transport covariance inflation.
#' @param bpe_not_applicable Behavior when subset diagnostics fail. One of `"na"`
#'   (default) to return NA estimates, or `"error"` to stop.
#' @param scale_instrument One of `"residual_sd"` (default) or `"none"`.
#' @param grid List controlling UCI bounds or other tuning parameters. For
#'   backward-compatible scalar UCI, if `grid$delta` is supplied and
#'   `grid$gmin`/`grid$gmax` are omitted, the package interprets `delta` as a
#'   direct-effect bound of `[-delta, +delta]` under the chosen
#'   `scale_instrument`. When `violation_pattern` is supplied, `grid$delta`
#'   instead refers to theta bounds over the pattern-scaled direct effect.
#' @param ... Reserved.
#'
#' @return Object of class `plausexog_fit`.
#'
#' @details
#' `spliv` implements patterned sensitivity analysis for exclusion violations in
#' IV designs with fixed effects or other residualization steps. Researchers can
#' supply a theoretically motivated `spliv_pattern()` object to specify where
#' direct effects of an instrument are expected to be larger or smaller, and the
#' package then scales LTZ/UCI sensitivity along that pattern.
#'
#' The package does not estimate an unrestricted direct-effect field. Instead,
#' researchers supply a structured pattern and ask whether conclusions survive
#' direct effects scaled along that pattern.
#'
#' In applied work, users should usually vary `delta` over a range with
#' `spliv_sensitivity_path()` rather than report one arbitrary sensitivity
#' value.
#'
#' Patterned sensitivity currently supports one endogenous treatment, one
#' excluded instrument, and one researcher-specified pattern at a time.
#'
#' Confirmatory BPE is not a subgroup-search procedure. The researcher must
#' supply a pre-specified `bpe_design()` object, the package validates that
#' subset, and BPE proceeds only if the confirmatory eligibility checks pass.
#'
#' The first-stage F-statistic is still reported for diagnostics, but
#' confirmatory BPE eligibility is determined by the pre-specification checks,
#' subset size, cluster count, residualized instrument variation, and a
#' first-stage equivalence interval.
#'
#' BPE reduced-form covariance is propagated as a full covariance matrix and can
#' optionally be inflated via `bpe_transport`.
#'
#' @name spliv
#' @examples
#' set.seed(1)
#' d <- data.frame(y = rnorm(60), x = rnorm(60), z = rnorm(60), w = rnorm(60))
#' fit <- spliv(y ~ x + w | z + w, d, method = "uci", delta = 0.1)
#' fit$estimates
#' @export
spliv <- function(formula,
                  data,
                  fe = NULL,
                  fe_engine = c("fixest", "lfe"),
                  vcov = c("iid", "hc1", "cluster"),
                  cluster = NULL,
                  method = c("ltz", "uci", "bpe"),
                  prior = NULL,
                  delta = NULL,
                  violation_pattern = NULL,
                  bpe = FALSE,
                  bpe_design = NULL,
                  bpe_spec = list(design = NULL, subset = NULL, z_names = NULL),
                  bpe_kappa = 1,
                  bpe_omega = NULL,
                  bpe_min_n_S = 2000,
                  bpe_min_clusters_S = 30,
                  bpe_max_F_S = NULL,
                  bpe_min_varZ_S = 1e-6,
                  bpe_equiv_margin = NULL,
                  bpe_equiv_level = 0.95,
                  bpe_transport = c("sampling", "none", "conservative"),
                  bpe_transport_kappa = 0,
                  bpe_not_applicable = c("na", "error"),
                  scale_instrument = c("residual_sd", "none"),
                  grid = list(),
                  ...) {
  .spliv_impl(
    formula = formula,
    data = data,
    fe = fe,
    fe_engine = fe_engine,
    vcov = vcov,
    cluster = cluster,
    method = method,
    prior = prior,
    delta = delta,
    violation_pattern = violation_pattern,
    bpe = bpe,
    bpe_design = bpe_design,
    bpe_spec = bpe_spec,
    bpe_kappa = bpe_kappa,
    bpe_omega = bpe_omega,
    bpe_min_n_S = bpe_min_n_S,
    bpe_min_clusters_S = bpe_min_clusters_S,
    bpe_max_F_S = bpe_max_F_S,
    bpe_min_varZ_S = bpe_min_varZ_S,
    bpe_equiv_margin = bpe_equiv_margin,
    bpe_equiv_level = bpe_equiv_level,
    bpe_transport = bpe_transport,
    bpe_transport_kappa = bpe_transport_kappa,
    bpe_not_applicable = bpe_not_applicable,
    scale_instrument = scale_instrument,
    grid = grid,
    ...
  )
}

.spliv_impl <- function(formula,
                        data,
                        fe = NULL,
                        fe_engine = c("fixest", "lfe"),
                        vcov = c("iid", "hc1", "cluster"),
                        cluster = NULL,
                        method = c("ltz", "uci", "bpe"),
                        prior = NULL,
                        delta = NULL,
                        violation_pattern = NULL,
                        bpe = FALSE,
                        bpe_design = NULL,
                        bpe_spec = list(design = NULL, subset = NULL, z_names = NULL),
                        bpe_kappa = 1,
                        bpe_omega = NULL,
                        bpe_min_n_S = 2000,
                        bpe_min_clusters_S = 30,
                        bpe_max_F_S = NULL,
                        bpe_min_varZ_S = 1e-6,
                        bpe_equiv_margin = NULL,
                        bpe_equiv_level = 0.95,
                        bpe_transport = c("sampling", "none", "conservative"),
                        bpe_transport_kappa = 0,
                        bpe_not_applicable = c("na", "error"),
                        scale_instrument = c("residual_sd", "none"),
                        grid = list(),
                        ...) {
  fe_engine <- match.arg(fe_engine)
  vcov <- match.arg(tolower(vcov), c("iid", "hc1", "cluster"))
  method_out <- match.arg(tolower(method), c("ltz", "uci", "bpe"))
  bpe_transport <- match.arg(tolower(bpe_transport), c("sampling", "none", "conservative"))
  bpe_not_applicable <- match.arg(tolower(bpe_not_applicable), c("na", "error"))
  scale_instrument <- match.arg(scale_instrument)

  method_core <- if (identical(method_out, "bpe")) "ltz" else method_out
  bpe_requested <- isTRUE(bpe) || identical(method_out, "bpe")

  if (bpe_requested && !identical(method_core, "ltz")) {
    stop("BPE is only available with LTZ (use method='bpe' or method='ltz', bpe=TRUE).")
  }

  extra_vars <- unique(c(all.vars(fe), all.vars(cluster)))
  parsed <- .iv_parse(formula, data, extra_vars = extra_vars)

  y <- parsed$y
  X <- parsed$X
  Z <- parsed$Z
  W <- .build_W_from_XZ(X, Z)
  data_cc <- data[parsed$keep, , drop = FALSE]

  if (!is.null(fe)) {
    X <- .drop_intercept_and_constants(X)
    Z <- .drop_intercept_and_constants(Z)
    W <- .drop_intercept_and_constants(W)
    if (ncol(X) == 0 || ncol(Z) == 0) {
      stop("After FE intercept/constant removal, X and Z must each have at least one column.")
    }

    if (fe_engine == "fixest") {
      dm <- demean_fixest(y = y, X = X, Z = Z, W = W, fe_fml = fe, data = data_cc)
    } else {
      fe_df <- .build_fe_frame(fe, data_cc)
      dm <- demean_lfe(y = y, X = X, Z = Z, W = W, fe_list = fe_df)
    }

    y <- dm$y
    X <- .drop_intercept_and_constants(dm$X)
    Z <- .drop_intercept_and_constants(dm$Z)
    W <- .drop_intercept_and_constants(dm$W)
    if (ncol(X) == 0 || ncol(Z) == 0) {
      stop("Demeaning removed all variation from X or Z.")
    }
  }

  cluster_id <- NULL
  if (vcov == "cluster") {
    if (is.null(cluster)) {
      stop("vcov='cluster' requires `cluster`.")
    }
    cluster_id <- .get_cluster_id(cluster, data_cc)
  }

  residualized_instrument_sd <- .instrument_residual_sds(
    X = X,
    Z = Z,
    W = W,
    fe_present = !is.null(fe)
  )

  bpe_diag <- NULL
  violation_pattern_state <- NULL
  violation_pattern_info <- NULL
  err_msgs <- character(0)
  warn_msgs <- character(0)

  if (!is.null(violation_pattern) && bpe_requested) {
    stop(
      "`violation_pattern` is not currently supported with confirmatory BPE. ",
      "Use patterned sensitivity with LTZ or UCI, or omit `violation_pattern` when `method = 'bpe'`."
    )
  }

  if (!is.null(violation_pattern)) {
    violation_pattern_state <- .prepare_violation_pattern(
      violation_pattern = violation_pattern,
      data = data,
      sample_idx = parsed$keep,
      X = X,
      Z = Z,
      residualized_instrument_sd = residualized_instrument_sd,
      scale_instrument = scale_instrument
    )
    violation_pattern_info <- violation_pattern_state$info
    warn_msgs <- c(warn_msgs, violation_pattern_state$warnings %||% character(0))
  }

  if (!is.null(bpe_omega)) {
    warn_msgs <- c(
      warn_msgs,
      "`bpe_omega` is deprecated; confirmatory BPE now uses the full reduced-form covariance and `bpe_transport`."
    )
  }
  if (!is.null(bpe_max_F_S)) {
    warn_msgs <- c(
      warn_msgs,
      "`bpe_max_F_S` is deprecated; the first-stage F-statistic is now reported for diagnostics only and does not determine confirmatory BPE eligibility."
    )
  }

  if (bpe_requested) {
    if (!is.null(prior)) {
      warn_msgs <- c(warn_msgs, "`prior` was supplied but BPE was requested; using the confirmatory BPE-learned prior.")
    }

    design_obj <- .coerce_bpe_design(
      bpe_design_arg = bpe_design,
      bpe_spec = bpe_spec,
      parsed = parsed,
      data = data
    )
    z_names <- bpe_spec$z_names %||% .default_uci_inst(X, Z)

    validation <- bpe_validate_design(
      formula = formula,
      data = data,
      design = design_obj,
      fe = fe,
      fe_engine = fe_engine,
      vcov = vcov,
      cluster = cluster,
      z_names = z_names,
      bpe_min_n_S = bpe_min_n_S,
      bpe_min_clusters_S = bpe_min_clusters_S,
      bpe_min_varZ_S = bpe_min_varZ_S,
      bpe_equiv_margin = bpe_equiv_margin,
      bpe_equiv_level = bpe_equiv_level,
      bpe_transport = bpe_transport,
      bpe_transport_kappa = bpe_transport_kappa,
      bpe_kappa = bpe_kappa,
      scale_instrument = scale_instrument
    )

    bpe_diag <- validation
    bpe_diag$not_applicable <- !isTRUE(validation$eligibility_passed)
    bpe_diag$subset_size <- validation$n_S
    bpe_diag$F_S <- as.numeric(validation$first_stage_f_statistic)
    bpe_diag$f_stat_S <- as.numeric(validation$first_stage_f_statistic)
    bpe_diag$gamma_S <- validation$reduced_form_direct_effect
    bpe_diag$SE_gamma_S <- sqrt(pmax(0, diag(validation$reduced_form_sampling_cov)))
    names(bpe_diag$SE_gamma_S) <- names(validation$reduced_form_direct_effect)
    bpe_diag$z_names <- validation$instrument

    warn_msgs <- c(warn_msgs, validation$warnings %||% character(0))
    bpe_prior <- .align_bpe_prior_to_Z(validation, Z)

    if (!isTRUE(validation$eligibility_passed)) {
      if (bpe_not_applicable == "error") {
        stop(validation$message)
      }

      err_msgs <- c(err_msgs, validation$message)
      warn_msgs <- c(warn_msgs, validation$message)

      return(.make_not_applicable_fit(
        call = match.call(),
        formula = formula,
        fe = fe,
        fe_engine = fe_engine,
        vcov = vcov,
        cluster = cluster,
        method_out = method_out,
        method_core = method_core,
        coef_names = colnames(X),
        prior = bpe_prior,
        grid = grid,
        bpe_diag = bpe_diag,
        msg = validation$message,
        y = y,
        X = X,
        Z = Z,
        W = W,
        cluster_id = cluster_id,
        scale_instrument = scale_instrument,
        residualized_instrument_sd = residualized_instrument_sd
      ))
    }

    prior <- bpe_prior
  }

  if (method_core == "ltz") {
    if (!is.null(violation_pattern_state)) {
      theta_prior <- .coerce_pattern_ltz_prior(prior = prior, delta = delta)
      ltz <- sp_ltz_mats(
        y = y,
        X = X,
        Z = Z,
        mu = theta_prior$mu,
        Omega = theta_prior$Omega,
        level = grid$level %||% 0.95,
        vcov = vcov,
        cluster_id = cluster_id,
        coef_names = colnames(X),
        inst_names = colnames(Z),
        direct_effect = violation_pattern_state$direct_effect,
        direct_effect_names = violation_pattern_state$direct_effect_names
      )

      fit <- list(
        call = match.call(),
        formula = formula,
        fe = fe,
        fe_engine = if (is.null(fe)) NULL else fe_engine,
        vcov = vcov,
        cluster = cluster,
        method = method_out,
        method_core = method_core,
        estimates = ltz$table,
        beta_hat = ltz$beta,
        beta_iv = ltz$beta_iv,
        vcov_beta = ltz$vcov,
        prior = list(mu = theta_prior$mu, Omega = theta_prior$Omega),
        mu_used = ltz$mu_used,
        Omega_used = as.matrix(theta_prior$Omega),
        A = ltz$A,
        A_mu = ltz$A_mu,
        beta_shift = ltz$beta_shift,
        grid = c(grid, list(delta = delta)),
        bpe = bpe_requested,
        violation_pattern = violation_pattern_info,
        diag = list(
          bpe = bpe_diag,
          violation_pattern = violation_pattern_info,
          ltz = list(
            beta_iv = ltz$beta_iv,
            mu_used = ltz$mu_used,
            A = ltz$A,
            A_mu = ltz$A_mu,
            beta_shift = ltz$beta_shift,
            direct_effect_names = ltz$direct_effect_names
          )
        ),
        bpe_diagnostics = bpe_diag,
        scale_instrument = scale_instrument,
        residualized_instrument_sd = residualized_instrument_sd,
        errors = err_msgs,
        warnings = unique(warn_msgs),
        internals = list(y = y, X = X, Z = Z, W = W, cluster_id = cluster_id)
      )
      class(fit) <- "plausexog_fit"
      return(fit)
    }

    if (!is.null(prior) && !is.null(delta)) {
      stop("Supply either `prior` or `delta` for LTZ, not both.")
    }
    if (is.null(prior) && !is.null(delta)) {
      prior <- .scalar_ltz_prior_from_delta(
        delta = delta,
        X = X,
        Z = Z,
        residualized_instrument_sd = residualized_instrument_sd,
        scale_instrument = scale_instrument
      )
    }
    if (is.null(prior)) {
      stop("LTZ requires `prior`, or request BPE via method='bpe' / bpe=TRUE with a confirmatory `bpe_design()`.")
    }

    mu <- prior$mu %||% prior$mu_hat
    Omega <- prior$Omega %||% prior$omega %||% prior$omega_hat
    if (is.null(mu) || is.null(Omega)) {
      stop("`prior` must include mu and Omega (or omega).")
    }

    ltz <- sp_ltz_mats(
      y = y,
      X = X,
      Z = Z,
      mu = mu,
      Omega = Omega,
      level = grid$level %||% 0.95,
      vcov = vcov,
      cluster_id = cluster_id,
      coef_names = colnames(X),
      inst_names = colnames(Z)
    )

    fit <- list(
      call = match.call(),
      formula = formula,
      fe = fe,
      fe_engine = if (is.null(fe)) NULL else fe_engine,
      vcov = vcov,
      cluster = cluster,
      method = method_out,
      method_core = method_core,
      estimates = ltz$table,
      beta_hat = ltz$beta,
      beta_iv = ltz$beta_iv,
      vcov_beta = ltz$vcov,
      prior = list(mu = mu, Omega = Omega),
      mu_used = ltz$mu_used,
      Omega_used = as.matrix(Omega),
      A = ltz$A,
      A_mu = ltz$A_mu,
      beta_shift = ltz$beta_shift,
      grid = c(grid, list(delta = delta)),
      bpe = bpe_requested,
      violation_pattern = violation_pattern_info,
      diag = list(
        bpe = bpe_diag,
        violation_pattern = violation_pattern_info,
        ltz = list(
          beta_iv = ltz$beta_iv,
          mu_used = ltz$mu_used,
          A = ltz$A,
          A_mu = ltz$A_mu,
          beta_shift = ltz$beta_shift
        )
      ),
      bpe_diagnostics = bpe_diag,
      scale_instrument = scale_instrument,
      residualized_instrument_sd = residualized_instrument_sd,
      errors = err_msgs,
      warnings = unique(warn_msgs),
      internals = list(y = y, X = X, Z = Z, W = W, cluster_id = cluster_id)
    )
    class(fit) <- "plausexog_fit"
    return(fit)
  }

  if (!is.null(delta) && (!is.null(grid$delta) || !is.null(grid$gmin) || !is.null(grid$gmax))) {
    stop("Supply either top-level `delta` or explicit `grid$delta`/`grid$gmin`/`grid$gmax`, not both.")
  }

  if (!is.null(violation_pattern_state)) {
    theta_delta <- delta %||% grid$delta %||% NULL
    steps <- grid$steps %||% grid$grid %||% 21
    if (length(steps) != 1) {
      stop("Patterned UCI currently supports a single grid length for one pattern coefficient.")
    }

    if (!is.null(theta_delta) && (is.null(grid$gmin) || is.null(grid$gmax))) {
      theta_delta <- as.numeric(theta_delta)
      if (length(theta_delta) != 1 || !is.finite(theta_delta) || theta_delta < 0) {
        stop("`delta` must be a non-negative finite scalar for patterned UCI.")
      }
      gmin <- grid$gmin %||% -theta_delta
      gmax <- grid$gmax %||% theta_delta
    } else {
      gmin <- grid$gmin %||% -1
      gmax <- grid$gmax %||% 1
    }
    if (length(gmin) != 1 || length(gmax) != 1 || !is.finite(gmin) || !is.finite(gmax)) {
      stop("Patterned UCI currently requires scalar theta bounds.")
    }

    gamma_grid <- .make_gamma_grid(gmin = as.numeric(gmin), gmax = as.numeric(gmax), steps = as.integer(steps))
    uci <- sp_uci_mats(
      y = y,
      X = X,
      Z = Z,
      inst_idx = 1L,
      gamma_grid = gamma_grid,
      level = grid$level %||% 0.95,
      vcov = vcov,
      cluster_id = cluster_id,
      coef_names = colnames(X),
      direct_effect = violation_pattern_state$direct_effect,
      direct_effect_names = violation_pattern_state$direct_effect_names
    )

    fit <- list(
      call = match.call(),
      formula = formula,
      fe = fe,
      fe_engine = if (is.null(fe)) NULL else fe_engine,
      vcov = vcov,
      cluster = cluster,
      method = method_out,
      method_core = method_core,
      estimates = uci,
      beta_hat = NA_real_,
      vcov_beta = NULL,
      prior = prior,
      mu_used = NULL,
      Omega_used = NULL,
      grid = c(
        grid,
        list(
          inst = violation_pattern_state$instrument,
          delta = theta_delta,
          gmin = as.numeric(gmin),
          gmax = as.numeric(gmax),
          steps = as.integer(steps),
          parameter = "theta"
        )
      ),
      bpe = bpe_requested,
      violation_pattern = violation_pattern_info,
      diag = list(
        bpe = bpe_diag,
        violation_pattern = violation_pattern_info
      ),
      bpe_diagnostics = bpe_diag,
      scale_instrument = scale_instrument,
      residualized_instrument_sd = residualized_instrument_sd,
      errors = err_msgs,
      warnings = unique(warn_msgs),
      internals = list(y = y, X = X, Z = Z, W = W, cluster_id = cluster_id)
    )
    class(fit) <- "plausexog_fit"
    return(fit)
  }

  if (!is.null(delta)) {
    grid$delta <- delta
  }

  inst <- grid$inst %||% .default_uci_inst(X, Z)
  if (is.character(inst)) {
    if (!all(inst %in% colnames(Z))) {
      stop("`grid$inst` contains names not in Z.")
    }
    inst_idx <- match(inst, colnames(Z))
  } else {
    inst_idx <- as.integer(inst)
    inst <- colnames(Z)[inst_idx]
  }

  p <- length(inst_idx)
  steps <- grid$steps %||% grid$grid %||% rep(21, p)
  if (length(steps) == 1) {
    steps <- rep(steps, p)
  }

  if (!is.null(grid$delta) && (is.null(grid$gmin) || is.null(grid$gmax))) {
    bounds <- .delta_bounds_from_scale(
      delta = grid$delta,
      inst = inst,
      residualized_instrument_sd = residualized_instrument_sd,
      scale_instrument = scale_instrument
    )
    gmin <- grid$gmin %||% bounds$gmin
    gmax <- grid$gmax %||% bounds$gmax
  } else {
    gmin <- grid$gmin %||% rep(-1, p)
    gmax <- grid$gmax %||% rep(1, p)
  }
  if (length(gmin) == 1) gmin <- rep(gmin, p)
  if (length(gmax) == 1) gmax <- rep(gmax, p)

  gamma_grid <- .make_gamma_grid(gmin = gmin, gmax = gmax, steps = steps)
  uci <- sp_uci_mats(
    y = y,
    X = X,
    Z = Z,
    inst_idx = inst_idx,
    gamma_grid = gamma_grid,
    level = grid$level %||% 0.95,
    vcov = vcov,
    cluster_id = cluster_id,
    coef_names = colnames(X)
  )

  fit <- list(
    call = match.call(),
    formula = formula,
    fe = fe,
    fe_engine = if (is.null(fe)) NULL else fe_engine,
    vcov = vcov,
    cluster = cluster,
    method = method_out,
    method_core = method_core,
    estimates = uci,
    beta_hat = NA_real_,
    vcov_beta = NULL,
    prior = prior,
    mu_used = NULL,
    Omega_used = NULL,
    grid = c(grid, list(inst = inst, gmin = gmin, gmax = gmax, steps = steps, parameter = "gamma")),
    bpe = bpe_requested,
    violation_pattern = violation_pattern_info,
    diag = list(
      bpe = bpe_diag,
      violation_pattern = violation_pattern_info
    ),
    bpe_diagnostics = bpe_diag,
    scale_instrument = scale_instrument,
    residualized_instrument_sd = residualized_instrument_sd,
    errors = err_msgs,
    warnings = unique(warn_msgs),
    internals = list(y = y, X = X, Z = Z, W = W, cluster_id = cluster_id)
  )
  class(fit) <- "plausexog_fit"
  fit
}

#' @export
print.plausexog_fit <- function(x, ...) {
  cat("spliv fit\n")
  cat("  method:", x$method, "\n")
  if (!is.null(x$fe)) {
    cat("  FE:", deparse(x$fe), "(engine:", x$fe_engine, ")\n")
  }
  if (!is.null(x$violation_pattern)) {
    cat("  violation pattern:", x$violation_pattern$name, "\n")
  }
  if (length(x$errors) > 0) {
    cat("  errors:", paste(x$errors, collapse = " | "), "\n")
  }
  print(x$estimates)
  invisible(x)
}
