#' Create a Confirmatory BPE Design Object
#'
#' Creates a researcher-specified design object for confirmatory BPE. The subset
#' must be justified outside the outcome analysis and should represent a
#' theoretically motivated instrument-inactive subset.
#'
#' @param name Short design name.
#' @param subset Subset specification. Accepts a one-sided formula, a
#'   `function(data)`, a logical vector of length `nrow(data)`, or a character
#'   string naming a logical column in `data`.
#' @param rationale Non-empty substantive justification for the subset.
#' @param variables_used Optional character vector describing the design
#'   variables used to define the subset. This is required for function-based
#'   subsets unless the package can safely infer the variables.
#' @param subset_type Optional short label such as `"theory_defined"` or
#'   `"design_based"`.
#' @param pre_specified Logical indicator for confirmatory use. Confirmatory BPE
#'   requires `TRUE`.
#' @param transportability_rationale Optional description of why the direct
#'   effect learned in the subset may be informative for the target sample.
#' @param notes Optional free-form notes.
#'
#' @return An object of class `"spliv_bpe_design"`.
#' @examples
#' design <- bpe_design(
#'   "Inactive subset", ~ inactive,
#'   rationale = "The treatment channel is absent in this subset.",
#'   variables_used = "inactive", pre_specified = TRUE
#' )
#' bpe_eval_subset(design, data.frame(inactive = c(TRUE, FALSE)))
#' @export
bpe_design <- function(name,
                       subset,
                       rationale,
                       variables_used = NULL,
                       subset_type = NULL,
                       pre_specified = TRUE,
                       transportability_rationale = NULL,
                       notes = NULL) {
  if (!is.character(name) || length(name) != 1 || !nzchar(trimws(name))) {
    stop("`name` must be a non-empty character scalar.")
  }
  if (missing(subset)) {
    stop("`subset` must be supplied.")
  }
  if (!is.character(rationale) || length(rationale) != 1) {
    stop("`rationale` must be a character scalar.")
  }
  if (!is.logical(pre_specified) || length(pre_specified) != 1 || is.na(pre_specified)) {
    stop("`pre_specified` must be TRUE or FALSE.")
  }

  subset_kind <- .bpe_subset_kind(subset)
  if (is.null(variables_used)) {
    variables_used <- .bpe_infer_variables(subset, subset_kind = subset_kind)
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
    subset = subset,
    rationale = rationale,
    variables_used = variables_used,
    subset_type = subset_type,
    pre_specified = pre_specified,
    transportability_rationale = transportability_rationale,
    notes = notes,
    subset_kind = subset_kind,
    created_at = Sys.time(),
    call = match.call()
  )
  class(out) <- "spliv_bpe_design"
  out
}

.bpe_subset_kind <- function(subset) {
  if (inherits(subset, "formula")) {
    if (length(subset) != 2) {
      stop("`subset` formulas for `bpe_design()` must be one-sided, e.g. `~ inactive_region`.")
    }
    return("formula")
  }
  if (is.function(subset)) {
    return("function")
  }
  if (is.logical(subset)) {
    return("logical")
  }
  if (is.character(subset) && length(subset) == 1 && nzchar(trimws(subset))) {
    return("column")
  }
  stop(
    "`subset` must be a one-sided formula, a function(data), a logical vector, ",
    "or a character string naming a logical column."
  )
}

.bpe_infer_variables <- function(subset, subset_kind = .bpe_subset_kind(subset)) {
  if (identical(subset_kind, "formula")) {
    return(all.vars(subset))
  }
  if (identical(subset_kind, "column")) {
    return(as.character(subset))
  }
  NULL
}

.bpe_design_variables <- function(design) {
  vars <- design$variables_used %||% .bpe_infer_variables(design$subset, subset_kind = design$subset_kind)
  if (is.null(vars)) {
    return(character(0))
  }
  unique(as.character(vars))
}

.bpe_eval_formula_subset <- function(formula, data) {
  expr <- formula[[2L]]
  eval(expr, envir = data, enclos = parent.frame())
}

.bpe_eval_function_subset <- function(fun, data) {
  fun(data)
}

.bpe_eval_column_subset <- function(col_name, data) {
  if (!col_name %in% names(data)) {
    stop("Subset column `", col_name, "` was not found in `data`.")
  }
  data[[col_name]]
}

#' Evaluate a Confirmatory BPE Subset
#'
#' Evaluates a `bpe_design()` object on a data frame and returns the resulting
#' logical indicator.
#'
#' @param design A `bpe_design()` object.
#' @param data Data frame used for evaluation.
#' @param max_na_share Maximum allowable share of `NA` values in the subset
#'   indicator before evaluation fails. The default is `0.05`.
#'
#' @return A logical vector of length `nrow(data)`.
#' @examples
#' design <- bpe_design("Inactive", ~ inactive,
#'   rationale = "The treatment channel is absent.")
#' bpe_eval_subset(design, data.frame(inactive = c(TRUE, FALSE)))
#' @export
bpe_eval_subset <- function(design, data, max_na_share = 0.05) {
  if (!inherits(design, "spliv_bpe_design")) {
    stop("`design` must inherit from class `spliv_bpe_design`.")
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  }
  if (!is.numeric(max_na_share) || length(max_na_share) != 1 || !is.finite(max_na_share) ||
      max_na_share < 0 || max_na_share > 1) {
    stop("`max_na_share` must be a number between 0 and 1.")
  }

  raw <- tryCatch(
    {
      if (identical(design$subset_kind, "formula")) {
        .bpe_eval_formula_subset(design$subset, data)
      } else if (identical(design$subset_kind, "function")) {
        .bpe_eval_function_subset(design$subset, data)
      } else if (identical(design$subset_kind, "logical")) {
        design$subset
      } else if (identical(design$subset_kind, "column")) {
        .bpe_eval_column_subset(design$subset, data)
      } else {
        stop("Unsupported subset kind: ", design$subset_kind)
      }
    },
    error = function(e) {
      stop("Failed to evaluate BPE subset `", design$name, "`: ", conditionMessage(e), call. = FALSE)
    }
  )

  if (!is.logical(raw)) {
    stop("BPE subset `", design$name, "` must evaluate to a logical vector.")
  }
  if (length(raw) != nrow(data)) {
    stop(
      "BPE subset `", design$name, "` returned length ", length(raw),
      " but `data` has ", nrow(data), " rows."
    )
  }

  na_idx <- is.na(raw)
  na_share <- mean(na_idx)
  if (all(na_idx)) {
    stop("BPE subset `", design$name, "` evaluated to all `NA`.")
  }
  if (na_share > max_na_share) {
    stop(
      "BPE subset `", design$name, "` produced too many `NA` values (",
      sprintf("%.1f%%", 100 * na_share), ")."
    )
  }

  warnings <- character(0)
  if (any(na_idx)) {
    warnings <- c(
      warnings,
      paste0(
        "Subset `", design$name, "` produced ", sum(na_idx),
        " `NA` value(s); those rows were treated as FALSE during BPE evaluation."
      )
    )
    raw[na_idx] <- FALSE
  }

  if (all(raw)) {
    stop("BPE subset `", design$name, "` selected all observations.")
  }
  if (!any(raw)) {
    stop("BPE subset `", design$name, "` selected no observations.")
  }

  structure(as.logical(raw), bpe_warnings = warnings, bpe_na_share = na_share)
}

.bpe_require_confirmatory_design <- function(design) {
  if (!inherits(design, "spliv_bpe_design")) {
    stop("Confirmatory BPE requires a `bpe_design()` object.")
  }
  if (!isTRUE(design$pre_specified)) {
    stop("Confirmatory BPE requires `pre_specified = TRUE` in `bpe_design()`.")
  }
  if (!is.character(design$rationale) || !nzchar(trimws(design$rationale))) {
    stop("Confirmatory BPE requires a non-empty `rationale` in `bpe_design()`.")
  }
  if (identical(design$subset_kind, "function") && !length(.bpe_design_variables(design))) {
    stop(
      "Function-based BPE designs must provide `variables_used` so the subset ",
      "definition can be audited for confirmatory use."
    )
  }
  invisible(design)
}

.bpe_build_regression_design <- function(main_var, z_name, Z_sub, W_sub, fe_present) {
  Z_interest <- as.matrix(Z_sub[, z_name, drop = FALSE])
  colnames(Z_interest) <- z_name

  W_sub <- if (is.null(W_sub)) {
    matrix(0, nrow = nrow(Z_interest), ncol = 0)
  } else {
    as.matrix(W_sub)
  }
  keep_w <- setdiff(colnames(W_sub), z_name)
  W_keep <- if (length(keep_w) == 0) {
    matrix(0, nrow = nrow(Z_interest), ncol = 0)
  } else {
    W_sub[, keep_w, drop = FALSE]
  }

  X_full <- cbind(Z_interest, W_keep)
  if (!fe_present) {
    X_full <- cbind("(Intercept)" = 1, X_full)
  }

  X_restricted <- W_keep
  if (!fe_present) {
    X_restricted <- cbind("(Intercept)" = 1, X_restricted)
  }

  list(
    y = as.numeric(main_var),
    X_full = as.matrix(X_full),
    X_restricted = as.matrix(X_restricted),
    z_name = z_name
  )
}

.bpe_residualize_vector <- function(vec, W_sub, fe_present) {
  vec <- as.numeric(vec)
  W_sub <- if (is.null(W_sub)) {
    matrix(0, nrow = length(vec), ncol = 0)
  } else {
    as.matrix(W_sub)
  }

  if (!ncol(W_sub)) {
    if (fe_present) {
      return(vec)
    }
    return(vec - mean(vec))
  }

  X <- W_sub
  if (!fe_present) {
    X <- cbind("(Intercept)" = 1, X)
  }
  fit <- .ols_fit_fast(vec, X)
  fit$resid
}

.bpe_residualized_sd <- function(Z_sub, W_sub, z_name, fe_present) {
  z_resid <- .bpe_residualize_vector(Z_sub[, z_name, drop = TRUE], W_sub = W_sub, fe_present = fe_present)
  stats::sd(as.numeric(z_resid))
}

.bpe_residualized_sd_from_vector <- function(vec, W_sub, fe_present) {
  vec_resid <- .bpe_residualize_vector(vec, W_sub = W_sub, fe_present = fe_present)
  stats::sd(as.numeric(vec_resid))
}

.bpe_confirmatory_first_stage_target <- function(X, Z) {
  x_names <- setdiff(colnames(X), "(Intercept)")
  z_names <- setdiff(colnames(Z), "(Intercept)")
  target_names <- setdiff(x_names, z_names)

  if (length(target_names) == 0) {
    stop(
      "Confirmatory BPE requires an identifiable endogenous treatment in `X` that is not also listed as an instrument."
    )
  }
  if (length(target_names) > 1) {
    stop(
      "Confirmatory BPE currently supports one endogenous treatment. ",
      "Please fit a model with exactly one endogenous regressor for confirmatory BPE."
    )
  }

  target_names[[1]]
}

.bpe_first_stage_f_statistic <- function(main_var, X_full, X_restricted) {
  fit_u <- .ols_fit_fast(main_var, X_full)
  rss_u <- sum(fit_u$resid^2)

  if (ncol(X_restricted) == 0) {
    rss_r <- sum(main_var^2)
  } else {
    fit_r <- .ols_fit_fast(main_var, X_restricted)
    rss_r <- sum(fit_r$resid^2)
  }

  q <- ncol(X_full) - ncol(X_restricted)
  n <- length(main_var)
  p_u <- ncol(X_full)
  den_df <- n - p_u
  if (q <= 0 || den_df <= 0 || !is.finite(rss_u) || rss_u <= 0) {
    return(NA_real_)
  }
  ((rss_r - rss_u) / q) / (rss_u / den_df)
}

.bpe_first_stage_diagnostics <- function(main_var,
                                         Z_sub,
                                         W_sub,
                                         z_name,
                                         fe_present,
                                         vcov,
                                         cluster_id = NULL,
                                         level = 0.95) {
  reg <- .bpe_build_regression_design(
    main_var = main_var,
    z_name = z_name,
    Z_sub = Z_sub,
    W_sub = W_sub,
    fe_present = fe_present
  )

  fit <- .iv_2sls_mats(
    y = reg$y,
    X = reg$X_full,
    Z = reg$X_full,
    vcov = vcov,
    cluster_id = cluster_id
  )

  z_pos <- match(z_name, colnames(reg$X_full))
  if (is.na(z_pos)) {
    stop("First-stage instrument `", z_name, "` was not found by column name.")
  }
  coef_hat <- as.numeric(fit$beta[z_pos])
  se_hat <- sqrt(pmax(0, diag(fit$vcov)[z_pos]))
  zcrit <- stats::qnorm(1 - (1 - level) / 2)
  ci <- c(coef_hat - zcrit * se_hat, coef_hat + zcrit * se_hat)
  f_stat <- .bpe_first_stage_f_statistic(
    main_var = reg$y,
    X_full = reg$X_full,
    X_restricted = reg$X_restricted
  )

  list(
    coefficient = setNames(coef_hat, z_name),
    se = setNames(se_hat, z_name),
    ci = matrix(ci, nrow = 1, dimnames = list(z_name, c("lower", "upper"))),
    f_statistic = setNames(f_stat, z_name),
    f_type = if (identical(vcov, "cluster")) "conventional_ols_diagnostic" else "conventional_ols_diagnostic"
  )
}

.first_stage_diag <- function(x_vec, Z_sub, W_sub, z_names, fe_present) {
  if (length(z_names) != 1) {
    stop("Legacy first-stage diagnostics currently support exactly one instrument.")
  }
  diag_obj <- .bpe_first_stage_diagnostics(
    main_var = x_vec,
    Z_sub = Z_sub,
    W_sub = W_sub,
    z_name = z_names[[1]],
    fe_present = fe_present,
    vcov = "iid",
    cluster_id = NULL,
    level = 0.95
  )
  list(
    var_z = stats::setNames(stats::var(as.numeric(Z_sub[, z_names[[1]], drop = TRUE])), z_names),
    pi_hat = diag_obj$coefficient,
    f_stat = diag_obj$f_statistic,
    q = length(z_names),
    n = length(x_vec)
  )
}

.bpe_normalize_equiv_margin <- function(bpe_equiv_margin, z_name) {
  if (missing(bpe_equiv_margin) || is.null(bpe_equiv_margin)) {
    stop("Confirmatory BPE requires a researcher-specified `bpe_equiv_margin`.")
  }
  if (is.numeric(bpe_equiv_margin) && length(bpe_equiv_margin) == 1) {
    if (!is.finite(bpe_equiv_margin) || bpe_equiv_margin <= 0) {
      stop("`bpe_equiv_margin` must be a positive finite scalar or a named numeric vector.")
    }
    return(setNames(as.numeric(bpe_equiv_margin), z_name))
  }
  if (is.numeric(bpe_equiv_margin) && !is.null(names(bpe_equiv_margin))) {
    if (!z_name %in% names(bpe_equiv_margin)) {
      stop("Named `bpe_equiv_margin` does not contain an entry for instrument `", z_name, "`.")
    }
    margin <- as.numeric(bpe_equiv_margin[[z_name]])
    if (!is.finite(margin) || margin <= 0) {
      stop("`bpe_equiv_margin` for instrument `", z_name, "` must be positive and finite.")
    }
    return(setNames(margin, z_name))
  }
  stop("`bpe_equiv_margin` must be a positive scalar or a named numeric vector.")
}

.bpe_transport_covariance <- function(vcov_mat,
                                      bpe_transport = c("none", "sampling", "conservative"),
                                      bpe_transport_kappa = 0) {
  bpe_transport <- match.arg(tolower(bpe_transport), c("none", "sampling", "conservative"))
  if (!is.matrix(vcov_mat)) {
    stop("`vcov_mat` must be a matrix.")
  }
  if (!is.numeric(bpe_transport_kappa) || length(bpe_transport_kappa) != 1 ||
      !is.finite(bpe_transport_kappa) || bpe_transport_kappa < 0) {
    stop("`bpe_transport_kappa` must be a non-negative finite scalar.")
  }

  if (bpe_transport %in% c("none", "sampling")) {
    return(list(
      vcov = as.matrix(vcov_mat),
      mode = bpe_transport,
      inflation = 1
    ))
  }

  list(
    vcov = (1 + bpe_transport_kappa) * as.matrix(vcov_mat),
    mode = bpe_transport,
    inflation = 1 + bpe_transport_kappa
  )
}

.bpe_collect_design_warnings <- function(eval_idx, audit_warnings) {
  unique(c(attr(eval_idx, "bpe_warnings", exact = TRUE) %||% character(0), audit_warnings))
}

.bpe_design_audit <- function(design, parsed, z_name) {
  vars_used <- .bpe_design_variables(design)
  x_names <- setdiff(colnames(parsed$X), "(Intercept)")
  z_names <- setdiff(colnames(parsed$Z), "(Intercept)")
  y_name <- all.vars(parsed$main_formula)[1]

  expr_text <- paste(deparse(design$subset), collapse = " ")
  lower_text <- tolower(expr_text)
  diagnostic_hit <- grepl("fitted|resid|residual|first_stage|f_stat|diagnostic|hatvalue|cooks", lower_text)

  warnings <- character(0)
  if (y_name %in% vars_used) {
    warnings <- c(
      warnings,
      paste0(
        "The BPE subset design references outcome variable `", y_name,
        "`. Confirmatory BPE subsets should not be defined using outcomes."
      )
    )
  }
  uses_endogenous <- intersect(vars_used, x_names)
  if (length(uses_endogenous)) {
    warnings <- c(
      warnings,
      paste0(
        "The BPE subset design references endogenous regressor(s): ",
        paste(uses_endogenous, collapse = ", "),
        ". Confirmatory designs should justify this choice carefully."
      )
    )
  }
  uses_instrument <- intersect(vars_used, z_names)
  if (length(uses_instrument)) {
    warnings <- c(
      warnings,
      paste0(
        "The BPE subset design references instrument variable(s): ",
        paste(uses_instrument, collapse = ", "),
        ". This audit warning does not by itself invalidate the design, but it should be justified."
      )
    )
  }
  if (diagnostic_hit) {
    warnings <- c(
      warnings,
      "The BPE subset definition appears to reference fitted-model diagnostics; confirmatory BPE should avoid data-driven subset construction."
    )
  }

  list(
    variables_used = vars_used,
    uses_outcome = y_name %in% vars_used,
    uses_endogenous = uses_endogenous,
    uses_instrument = uses_instrument,
    diagnostic_hit = diagnostic_hit,
    warnings = unique(warnings),
    expression = expr_text,
    instrument = z_name
  )
}

.bpe_full_sample_sd <- function(parsed, data, fe, fe_engine, z_name) {
  full_idx <- rep(TRUE, nrow(data))
  mats <- .bpe_build_subset_mats(
    parsed = parsed,
    data = data,
    subset_idx_full = full_idx,
    fe = fe,
    fe_engine = fe_engine
  )
  if (!isTRUE(mats$ok) || !z_name %in% colnames(mats$Z)) {
    return(NA_real_)
  }
  .bpe_residualized_sd(
    Z_sub = mats$Z,
    W_sub = mats$W,
    z_name = z_name,
    fe_present = !is.null(fe)
  )
}

.bpe_validation_from_parsed <- function(parsed,
                                        data,
                                        design,
                                        fe = NULL,
                                        fe_engine = c("fixest", "lfe"),
                                        vcov = c("iid", "hc1", "cluster"),
                                        cluster = NULL,
                                        z_names = NULL,
                                        bpe_min_n_S = 2000,
                                        bpe_min_clusters_S = 30,
                                        bpe_min_varZ_S = 1e-6,
                                        bpe_equiv_margin,
                                        bpe_equiv_level = 0.95,
                                        bpe_transport = c("none", "sampling", "conservative"),
                                        bpe_transport_kappa = 0,
                                        bpe_kappa = 1,
                                        scale_instrument = c("residual_sd", "none")) {
  .bpe_require_confirmatory_design(design)

  fe_engine <- match.arg(fe_engine)
  vcov <- match.arg(tolower(vcov), c("iid", "hc1", "cluster"))
  scale_instrument <- match.arg(scale_instrument)
  if (!is.numeric(bpe_min_n_S) || length(bpe_min_n_S) != 1 || !is.finite(bpe_min_n_S) || bpe_min_n_S <= 0) {
    stop("`bpe_min_n_S` must be a positive finite scalar.")
  }
  if (!is.numeric(bpe_min_clusters_S) || length(bpe_min_clusters_S) != 1 || !is.finite(bpe_min_clusters_S) || bpe_min_clusters_S < 2) {
    stop("`bpe_min_clusters_S` must be a finite scalar greater than or equal to 2.")
  }
  if (!is.numeric(bpe_min_varZ_S) || length(bpe_min_varZ_S) != 1 || !is.finite(bpe_min_varZ_S) || bpe_min_varZ_S <= 0) {
    stop("`bpe_min_varZ_S` must be a positive finite scalar.")
  }
  if (!is.numeric(bpe_equiv_level) || length(bpe_equiv_level) != 1 || !is.finite(bpe_equiv_level) ||
      bpe_equiv_level <= 0 || bpe_equiv_level >= 1) {
    stop("`bpe_equiv_level` must be a number strictly between 0 and 1.")
  }
  if (!is.numeric(bpe_kappa) || length(bpe_kappa) != 1 || !is.finite(bpe_kappa) || bpe_kappa <= 0) {
    stop("`bpe_kappa` must be a positive finite scalar.")
  }

  cluster_id_cc <- NULL
  if (vcov == "cluster") {
    data_cc <- data[parsed$keep, , drop = FALSE]
    cluster_id_cc <- .get_cluster_id(cluster, data_cc)
  }

  z_names <- z_names %||% .default_uci_inst(parsed$X, parsed$Z)
  z_names <- as.character(z_names)
  if (length(z_names) != 1) {
    stop("Confirmatory BPE currently supports exactly one instrument.")
  }
  z_name <- z_names[[1]]
  if (!z_name %in% colnames(parsed$Z)) {
    stop("Requested BPE instrument `", z_name, "` was not found in the parsed instrument matrix.")
  }
  fs_target <- .bpe_confirmatory_first_stage_target(parsed$X, parsed$Z)

  subset_idx_full <- bpe_eval_subset(design, data = data)
  audit <- .bpe_design_audit(design, parsed = parsed, z_name = z_name)
  warnings <- .bpe_collect_design_warnings(subset_idx_full, audit$warnings)

  mats_S <- .bpe_build_subset_mats(
    parsed = parsed,
    data = data,
    subset_idx_full = subset_idx_full,
    fe = fe,
    fe_engine = fe_engine
  )
  if (!isTRUE(mats_S$ok)) {
    out <- list(
      design_name = design$name,
      rationale = design$rationale,
      subset_type = design$subset_type,
      variables_used = .bpe_design_variables(design),
      pre_specified = design$pre_specified,
      transportability_rationale = design$transportability_rationale,
      notes = design$notes,
      n_S = 0,
      share_S = 0,
      G_S = if (identical(vcov, "cluster")) 0 else NULL,
      varZ_S = setNames(NA_real_, z_name),
      residualized_instrument_sd_S = setNames(NA_real_, z_name),
      residualized_instrument_sd = setNames(.bpe_full_sample_sd(parsed, data, fe, fe_engine, z_name), z_name),
      residualized_treatment_sd_S = setNames(NA_real_, fs_target),
      first_stage_coefficient = setNames(NA_real_, z_name),
      first_stage_se = setNames(NA_real_, z_name),
      first_stage_ci = matrix(NA_real_, nrow = 1, dimnames = list(z_name, c("lower", "upper"))),
      first_stage_f_statistic = setNames(NA_real_, z_name),
      first_stage_f_type = "conventional_ols_diagnostic",
      first_stage_effect_one_residual_sd_Z = setNames(NA_real_, z_name),
      standardized_first_stage_effect = setNames(NA_real_, z_name),
      equivalence_margin = .bpe_normalize_equiv_margin(bpe_equiv_margin, z_name),
      equivalence_level = bpe_equiv_level,
      equivalence_passed = FALSE,
      eligibility_passed = FALSE,
      eligibility_checks = list(
        pre_specified = TRUE,
        rationale = TRUE,
        minimum_n = FALSE,
        minimum_clusters = if (identical(vcov, "cluster")) FALSE else TRUE,
        residual_variation = FALSE,
        equivalence = FALSE
      ),
      reduced_form_direct_effect = setNames(NA_real_, z_name),
      reduced_form_direct_effect_cov = matrix(NA_real_, 1, 1, dimnames = list(z_name, z_name)),
      reduced_form_sampling_cov = matrix(NA_real_, 1, 1, dimnames = list(z_name, z_name)),
      transport_covariance = matrix(NA_real_, 1, 1, dimnames = list(z_name, z_name)),
      transport_mode = match.arg(tolower(bpe_transport), c("none", "sampling", "conservative")),
      transport_uncertainty_inflation = NA_real_,
      prior_mu_sub = setNames(NA_real_, z_name),
      prior_Omega_sub = matrix(NA_real_, 1, 1, dimnames = list(z_name, z_name)),
      prior_mu_full = rep(NA_real_, ncol(parsed$Z)),
      prior_Omega_full = matrix(NA_real_, ncol(parsed$Z), ncol(parsed$Z), dimnames = list(colnames(parsed$Z), colnames(parsed$Z))),
      subset_idx_full = subset_idx_full,
      design_audit = audit,
      warnings = unique(c(warnings, mats_S$reason)),
      message = mats_S$reason,
      instrument = z_name,
      scale_instrument = scale_instrument
    )
    class(out) <- "spliv_bpe_validation"
    return(out)
  }

  cluster_S <- if (identical(vcov, "cluster")) cluster_id_cc[mats_S$idx_cc] else NULL
  G_S <- if (is.null(cluster_S)) NULL else length(unique(as.character(cluster_S)))
  share_S <- nrow(mats_S$Z) / length(parsed$y)

  residual_sd_S <- .bpe_residualized_sd(
    Z_sub = mats_S$Z,
    W_sub = mats_S$W,
    z_name = z_name,
    fe_present = !is.null(fe)
  )
  varZ_S <- residual_sd_S^2
  residual_sd_full <- .bpe_full_sample_sd(parsed, data, fe, fe_engine, z_name)
  residual_sd_x_S <- .bpe_residualized_sd_from_vector(
    mats_S$X[, fs_target, drop = TRUE],
    W_sub = mats_S$W,
    fe_present = !is.null(fe)
  )

  fs <- .bpe_first_stage_diagnostics(
    main_var = mats_S$X[, fs_target, drop = TRUE],
    Z_sub = mats_S$Z,
    W_sub = mats_S$W,
    z_name = z_name,
    fe_present = !is.null(fe),
    vcov = vcov,
    cluster_id = cluster_S,
    level = bpe_equiv_level
  )
  margin <- .bpe_normalize_equiv_margin(bpe_equiv_margin, z_name)
  ci_row <- fs$ci[z_name, , drop = TRUE]
  equivalence_passed <- isTRUE(ci_row["lower"] >= -margin[[z_name]] && ci_row["upper"] <= margin[[z_name]])
  first_stage_effect_one_sd_z <- unname(fs$coefficient[[z_name]]) * residual_sd_S
  standardized_first_stage_effect <- if (is.finite(residual_sd_x_S) && residual_sd_x_S > 0) {
    first_stage_effect_one_sd_z / residual_sd_x_S
  } else {
    NA_real_
  }

  rf_S <- .bpe_reduced_form_gamma(
    y = mats_S$y,
    Z = mats_S$Z,
    W = mats_S$W,
    z_names = z_name,
    fe_present = !is.null(fe),
    vcov = vcov,
    cluster_id = cluster_S
  )
  transport <- .bpe_transport_covariance(
    vcov_mat = rf_S$vcov,
    bpe_transport = bpe_transport,
    bpe_transport_kappa = bpe_transport_kappa
  )
  prior_Omega_sub <- bpe_kappa * transport$vcov
  prior_full <- .embed_prior_by_names(
    inst_names = colnames(parsed$Z),
    z_names = z_name,
    mu_hat = rf_S$gamma,
    omega_hat = prior_Omega_sub
  )

  minimum_n_passed <- nrow(mats_S$Z) >= bpe_min_n_S
  minimum_clusters_passed <- if (is.null(G_S)) TRUE else G_S >= bpe_min_clusters_S
  residual_variation_passed <- is.finite(varZ_S) && varZ_S >= bpe_min_varZ_S
  eligibility_passed <- all(
    c(minimum_n_passed, minimum_clusters_passed, residual_variation_passed, equivalence_passed)
  )

  out <- list(
    design_name = design$name,
    rationale = design$rationale,
    subset_type = design$subset_type,
    variables_used = .bpe_design_variables(design),
    pre_specified = design$pre_specified,
    transportability_rationale = design$transportability_rationale,
    notes = design$notes,
    n_S = nrow(mats_S$Z),
    share_S = share_S,
    G_S = G_S,
    varZ_S = setNames(varZ_S, z_name),
    residualized_instrument_sd_S = setNames(residual_sd_S, z_name),
    residualized_instrument_sd = setNames(residual_sd_full, z_name),
    residualized_treatment_sd_S = setNames(residual_sd_x_S, fs_target),
    first_stage_coefficient = fs$coefficient,
    first_stage_se = fs$se,
    first_stage_ci = fs$ci,
    first_stage_f_statistic = fs$f_statistic,
    first_stage_f_type = fs$f_type,
    first_stage_effect_one_residual_sd_Z = setNames(first_stage_effect_one_sd_z, z_name),
    standardized_first_stage_effect = setNames(standardized_first_stage_effect, z_name),
    equivalence_margin = margin,
    equivalence_level = bpe_equiv_level,
    equivalence_passed = equivalence_passed,
    eligibility_passed = eligibility_passed,
    eligibility_checks = list(
      pre_specified = TRUE,
      rationale = TRUE,
      minimum_n = minimum_n_passed,
      minimum_clusters = minimum_clusters_passed,
      residual_variation = residual_variation_passed,
      equivalence = equivalence_passed
    ),
    reduced_form_direct_effect = setNames(as.numeric(rf_S$gamma), z_name),
    reduced_form_direct_effect_cov = rf_S$vcov,
    reduced_form_sampling_cov = rf_S$vcov,
    transport_covariance = transport$vcov,
    transport_mode = transport$mode,
    transport_uncertainty_inflation = transport$inflation,
    prior_mu_sub = setNames(as.numeric(rf_S$gamma), z_name),
    prior_Omega_sub = prior_Omega_sub,
    prior_mu_full = prior_full$mu,
    prior_Omega_full = prior_full$Omega,
    subset_idx_full = subset_idx_full,
    design_audit = audit,
    warnings = warnings,
    message = if (eligibility_passed) "" else .bpe_eligibility_message(
      n_S = nrow(mats_S$Z),
      G_S = G_S,
      varZ_S = varZ_S,
      margin = margin[[z_name]],
      ci_lower = ci_row["lower"],
      ci_upper = ci_row["upper"]
    ),
    first_stage_target = fs_target,
    instrument = z_name,
    scale_instrument = scale_instrument
  )
  class(out) <- "spliv_bpe_validation"
  out
}

.bpe_eligibility_message <- function(n_S, G_S, varZ_S, margin, ci_lower, ci_upper) {
  cluster_text <- if (is.null(G_S)) "not clustered" else paste0("G_S=", G_S)
  sprintf(
    paste0(
      "BPE eligibility failed for the proposed subset S ",
      "(n_S=%s, %s, varZ_S=%.6g, first-stage CI=[%.6g, %.6g], equivalence margin=[%.6g, %.6g])."
    ),
    as.character(round(n_S)),
    cluster_text,
    varZ_S,
    ci_lower,
    ci_upper,
    -margin,
    margin
  )
}

#' Validate a Confirmatory BPE Design
#'
#' Runs the confirmatory BPE eligibility diagnostics for a pre-specified design
#' without fitting the final SPLIV model.
#'
#' @param formula IV formula `y ~ X | Z`.
#' @param data Data frame.
#' @param design A `bpe_design()` object.
#' @param fe Optional one-sided formula of fixed effects.
#' @param fe_engine FE demeaning engine, one of `"fixest"` or `"lfe"`.
#' @param vcov One of `"iid"`, `"hc1"`, or `"cluster"`.
#' @param cluster Cluster ids or one-sided formula when `vcov = "cluster"`.
#' @param z_names Optional instrument name for BPE. Confirmatory BPE currently
#'   supports exactly one instrument.
#' @param bpe_min_n_S Minimum subset size threshold.
#' @param bpe_min_clusters_S Minimum number of clusters required in the subset
#'   when clustered covariance is used.
#' @param bpe_min_varZ_S Minimum residualized instrument variance required in
#'   the subset.
#' @param bpe_equiv_margin Researcher-specified equivalence margin for the
#'   first-stage coefficient. Eligibility is based on the first-stage
#'   equivalence interval, not on the first-stage F-statistic.
#' @param bpe_equiv_level Confidence level used for the first-stage equivalence
#'   interval.
#' @param bpe_transport One of `"none"`, `"sampling"`, or `"conservative"`.
#'   Transportability is an assumption reflected in the reported covariance; it
#'   is not established by the subset itself.
#' @param bpe_transport_kappa Non-negative scalar controlling the conservative
#'   transport covariance inflation.
#' @param bpe_kappa Positive scalar multiplier applied to the transported BPE
#'   covariance before it is embedded into the LTZ prior.
#' @param scale_instrument One of `"residual_sd"` or `"none"`.
#'
#' @return An object of class `"spliv_bpe_validation"` containing design
#'   metadata, subset diagnostics, first-stage equivalence diagnostics, reduced-
#'   form direct-effect estimates, and covariance components for confirmatory
#'   BPE.
#' @examples
#' set.seed(2)
#' d <- data.frame(
#'   y = rnorm(80), x = rnorm(80), z = rnorm(80),
#'   inactive = rep(c(TRUE, FALSE), each = 40)
#' )
#' design <- bpe_design("Inactive", ~ inactive,
#'   rationale = "The treatment channel is absent.")
#' bpe_validate_design(y ~ x | z, d, design,
#'   bpe_min_n_S = 20, bpe_equiv_margin = 1)
#' @export
bpe_validate_design <- function(formula,
                                data,
                                design,
                                fe = NULL,
                                fe_engine = c("fixest", "lfe"),
                                vcov = c("iid", "hc1", "cluster"),
                                cluster = NULL,
                                z_names = NULL,
                                bpe_min_n_S = 2000,
                                bpe_min_clusters_S = 30,
                                bpe_min_varZ_S = 1e-6,
                                bpe_equiv_margin,
                                bpe_equiv_level = 0.95,
                                bpe_transport = c("none", "sampling", "conservative"),
                                bpe_transport_kappa = 0,
                                bpe_kappa = 1,
                                scale_instrument = c("residual_sd", "none")) {
  extra_vars <- unique(c(all.vars(fe), all.vars(cluster)))
  parsed <- .iv_parse(formula, data, extra_vars = extra_vars)
  .bpe_validation_from_parsed(
    parsed = parsed,
    data = data,
    design = design,
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
}
