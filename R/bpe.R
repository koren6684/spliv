.ols_fit_fast <- function(y, X) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  XX <- crossprod(X)
  XXinv <- .safe_solve(XX)
  beta <- as.numeric(XXinv %*% crossprod(X, y))
  resid <- y - as.numeric(X %*% beta)
  n <- length(y)
  p <- ncol(X)
  s2 <- sum(resid^2) / max(1, n - p)
  V <- s2 * XXinv
  list(coef = beta, vcov = V, resid = resid, n = n, p = p)
}

.bpe_build_subset_mats <- function(parsed, data, subset_idx_full, fe = NULL, fe_engine = "fixest") {
  if (!is.logical(subset_idx_full) || length(subset_idx_full) != nrow(data)) {
    stop("`subset_idx_full` must be logical with length nrow(data).")
  }

  idx_cc <- as.logical(subset_idx_full[parsed$keep])
  if (!any(idx_cc)) {
    return(list(ok = FALSE, reason = "subset selects no complete-case observations."))
  }

  y_raw <- parsed$y
  X_raw <- parsed$X
  Z_raw <- parsed$Z
  W_raw <- .build_W_from_XZ(X_raw, Z_raw)
  data_cc <- data[parsed$keep, , drop = FALSE]

  y <- y_raw[idx_cc]
  X <- X_raw[idx_cc, , drop = FALSE]
  Z <- Z_raw[idx_cc, , drop = FALSE]
  W <- if (is.null(W_raw) || ncol(W_raw) == 0) W_raw else W_raw[idx_cc, , drop = FALSE]
  dsub <- data_cc[idx_cc, , drop = FALSE]

  if (!is.null(fe)) {
    X <- .drop_intercept_and_constants(X)
    Z <- .drop_intercept_and_constants(Z)
    W <- .drop_intercept_and_constants(W)

    if (ncol(X) == 0 || ncol(Z) == 0) {
      return(list(ok = FALSE, reason = "subset FE demeaning removed all variation from X or Z."))
    }

    if (fe_engine == "fixest") {
      dm <- .demean_fixest(y = y, X = X, Z = Z, W = W, fe_fml = fe, data = dsub)
    } else {
      fe_df <- .build_fe_frame(fe, dsub)
      dm <- .demean_lfe(y = y, X = X, Z = Z, W = W, fe_list = fe_df)
    }

    y <- dm$y
    X <- .drop_intercept_and_constants(dm$X)
    Z <- .drop_intercept_and_constants(dm$Z)
    W <- .drop_intercept_and_constants(dm$W)

    if (ncol(X) == 0 || ncol(Z) == 0) {
      return(list(ok = FALSE, reason = "subset FE demeaning removed all variation from X or Z."))
    }
  }

  list(
    ok = TRUE,
    y = as.numeric(y),
    X = as.matrix(X),
    Z = as.matrix(Z),
    W = if (is.null(W)) matrix(0, nrow = length(y), ncol = 0) else as.matrix(W),
    idx_cc = idx_cc
  )
}

.bpe_reduced_form_gamma <- function(y, Z, W, z_names, fe_present = FALSE, vcov = "iid", cluster_id = NULL) {
  if (!all(z_names %in% colnames(Z))) {
    stop("Some `z_names` are not in Z for reduced-form diagnostics.")
  }

  Z_interest <- Z[, z_names, drop = FALSE]
  keep_w <- setdiff(colnames(W), z_names)
  W_keep <- if (length(keep_w) == 0) matrix(0, nrow = nrow(Z), ncol = 0) else W[, keep_w, drop = FALSE]

  Xrf <- cbind(Z_interest, W_keep)
  if (!fe_present) {
    Xrf <- cbind("(Intercept)" = 1, Xrf)
  }

  fit <- .iv_2sls_mats(y, X = Xrf, Z = Xrf, vcov = vcov, cluster_id = cluster_id)
  z_pos <- match(z_names, colnames(Xrf))
  se <- sqrt(pmax(0, diag(fit$vcov)[z_pos]))
  if (any(!is.finite(se))) {
    return(list(
      gamma = rep(NA_real_, length(z_names)),
      se = rep(NA_real_, length(z_names)),
      vcov = matrix(NA_real_, length(z_names), length(z_names))
    ))
  }

  list(
    gamma = stats::setNames(as.numeric(fit$beta[z_pos]), z_names),
    se = stats::setNames(as.numeric(se), z_names),
    vcov = {
      vc <- fit$vcov[z_pos, z_pos, drop = FALSE]
      dimnames(vc) <- list(z_names, z_names)
      vc
    }
  )
}

.bpe_subset_diagnostics_from_parsed <- function(parsed,
                                                data,
                                                subset_idx_full,
                                                fe = NULL,
                                                fe_engine = "fixest",
                                                z_names = NULL,
                                                vcov = "iid",
                                                cluster_id_cc = NULL) {
  mats_S <- .bpe_build_subset_mats(
    parsed = parsed,
    data = data,
    subset_idx_full = subset_idx_full,
    fe = fe,
    fe_engine = fe_engine
  )

  if (!isTRUE(mats_S$ok)) {
    return(list(
      ok = FALSE,
      message = mats_S$reason,
      n_S = 0,
      varZ_S = NA_real_,
      F_S = NA_real_,
      pi_hat_S = NA_real_,
      gamma_S = NA_real_,
      SE_gamma_S = NA_real_,
      gamma_notS = NA_real_,
      SE_gamma_notS = NA_real_,
      z_names = character(0),
      vcov_used = vcov
    ))
  }

  z_names_eff <- z_names %||% .default_uci_inst(mats_S$X, mats_S$Z)
  if (!all(z_names_eff %in% colnames(mats_S$Z))) {
    return(list(
      ok = FALSE,
      message = "Requested BPE instrument(s) are not available in subset S after preprocessing.",
      n_S = nrow(mats_S$Z),
      varZ_S = NA_real_,
      F_S = NA_real_,
      pi_hat_S = NA_real_,
      gamma_S = NA_real_,
      SE_gamma_S = NA_real_,
      gamma_notS = NA_real_,
      SE_gamma_notS = NA_real_,
      z_names = z_names_eff,
      vcov_used = vcov
    ))
  }

  fs_target <- .first_stage_target(mats_S$X, mats_S$Z)
  fs <- .first_stage_diag(
    x_vec = mats_S$X[, fs_target, drop = TRUE],
    Z_sub = mats_S$Z,
    W_sub = mats_S$W,
    z_names = z_names_eff,
    fe_present = !is.null(fe)
  )

  cluster_S <- if (vcov == "cluster") {
    if (is.null(cluster_id_cc)) stop("vcov='cluster' but cluster_id_cc is NULL")
    cluster_id_cc[mats_S$idx_cc]
  } else {
    NULL
  }

  rf_S <- .bpe_reduced_form_gamma(
    y = mats_S$y,
    Z = mats_S$Z,
    W = mats_S$W,
    z_names = z_names_eff,
    fe_present = !is.null(fe),
    vcov = vcov,
    cluster_id = cluster_S
  )

  if (any(!is.finite(rf_S$se))) {
    return(list(
      ok = FALSE,
      message = "BPE reduced-form estimation produced non-finite standard errors.",
      n_S = nrow(mats_S$Z),
      varZ_S = as.numeric(fs$var_z),
      F_S = as.numeric(fs$f_stat),
      pi_hat_S = as.numeric(fs$pi_hat),
      gamma_S = NA_real_,
      SE_gamma_S = NA_real_,
      gamma_notS = NA_real_,
      SE_gamma_notS = NA_real_,
      z_names = z_names_eff,
      vcov_used = vcov
    ))
  }

  idx_notS <- !as.logical(subset_idx_full)
  mats_notS <- .bpe_build_subset_mats(
    parsed = parsed,
    data = data,
    subset_idx_full = idx_notS,
    fe = fe,
    fe_engine = fe_engine
  )

  if (isTRUE(mats_notS$ok) && all(z_names_eff %in% colnames(mats_notS$Z))) {
    cluster_notS <- if (vcov == "cluster") {
      cluster_id_cc[mats_notS$idx_cc]
    } else {
      NULL
    }

    rf_notS <- .bpe_reduced_form_gamma(
      y = mats_notS$y,
      Z = mats_notS$Z,
      W = mats_notS$W,
      z_names = z_names_eff,
      fe_present = !is.null(fe),
      vcov = vcov,
      cluster_id = cluster_notS
    )
    gamma_notS <- rf_notS$gamma
    se_notS <- rf_notS$se
    if (any(!is.finite(se_notS))) {
      gamma_notS <- rep(NA_real_, length(z_names_eff))
      se_notS <- rep(NA_real_, length(z_names_eff))
    }
  } else {
    gamma_notS <- rep(NA_real_, length(z_names_eff))
    se_notS <- rep(NA_real_, length(z_names_eff))
  }

  list(
    ok = TRUE,
    message = NULL,
    n_S = nrow(mats_S$Z),
    varZ_S = setNames(as.numeric(fs$var_z), z_names_eff),
    F_S = as.numeric(fs$f_stat),
    pi_hat_S = setNames(as.numeric(fs$pi_hat), z_names_eff),
    gamma_S = setNames(as.numeric(rf_S$gamma), z_names_eff),
    SE_gamma_S = setNames(as.numeric(rf_S$se), z_names_eff),
    gamma_notS = setNames(as.numeric(gamma_notS), z_names_eff),
    SE_gamma_notS = setNames(as.numeric(se_notS), z_names_eff),
    z_names = z_names_eff,
    mats_S = mats_S,
    mats_notS = mats_notS,
    vcov_used = vcov
  )
}

.bpe_exploration_warning_text <- function() {
  paste(
    "This function is exploratory. It should not be used for confirmatory BPE",
    "inference. Confirmatory BPE requires a pre-specified bpe_design() object",
    "with a substantive rationale."
  )
}

.bpe_normalize_exploratory_rules <- function(rules) {
  if (missing(rules) || is.null(rules) || (is.list(rules) && length(rules) == 0)) {
    stop(
      "`rules` must supply one or more candidate subset definitions as one-sided formulas, ",
      "functions, logical vectors, logical-column names, or named lists with `subset=`."
    )
  }

  if (inherits(rules, "formula") || is.function(rules) || is.logical(rules) ||
      (is.character(rules) && length(rules) == 1)) {
    rules <- list(candidate_1 = rules)
  } else if (is.character(rules) && length(rules) > 1) {
    rules <- stats::setNames(as.list(rules), rules)
  } else if (!is.list(rules)) {
    stop(
      "`rules` must be a formula, function, logical vector, character column name, ",
      "or a list of such candidate definitions."
    )
  }

  out <- vector("list", length(rules))
  nm <- names(rules)
  for (i in seq_along(rules)) {
    candidate_i <- rules[[i]]
    name_i <- nm[[i]] %||% paste0("candidate_", i)
    variables_used_i <- NULL

    if (is.list(candidate_i) &&
        !inherits(candidate_i, "formula") &&
        !is.function(candidate_i) &&
        !is.logical(candidate_i) &&
        !(is.character(candidate_i) && length(candidate_i) == 1)) {
      if (is.null(candidate_i$subset)) {
        stop("Exploratory rule lists must include a `subset` component.")
      }
      name_i <- candidate_i$name %||% name_i
      variables_used_i <- candidate_i$variables_used %||% NULL
      candidate_i <- candidate_i$subset
    }

    out[[i]] <- list(name = name_i, subset = candidate_i, variables_used = variables_used_i)
  }
  out
}

.bpe_exploratory_subset <- function(candidate, data) {
  design_i <- bpe_design(
    name = candidate$name,
    subset = candidate$subset,
    rationale = "Exploratory candidate subset. Not valid for confirmatory BPE.",
    variables_used = candidate$variables_used %||% NULL,
    subset_type = "exploratory_candidate",
    pre_specified = FALSE
  )
  idx <- bpe_eval_subset(design_i, data = data)
  list(
    design = design_i,
    subset = as.logical(idx),
    warnings = attr(idx, "bpe_warnings", exact = TRUE) %||% character(0)
  )
}

#' Explore Candidate BPE Subsets
#'
#' Explores user-supplied candidate subset definitions and reports diagnostics
#' for transparent exploratory work. This function is not valid confirmatory BPE
#' by itself. Confirmatory BPE requires a pre-specified `bpe_design()` object
#' supplied to `spliv()` or `bpe_validate_design()`.
#'
#' @param data Data frame used for estimation.
#' @param spec Either an IV formula (`y ~ X | Z`) or a list with at least
#'   `formula`, and optional `fe`, `fe_engine`, and `z_names`.
#' @param rules Candidate subset definitions. Each candidate may be a one-sided
#'   formula evaluated in `data`, a `function(data)` returning a logical vector,
#'   a character string naming a logical column in `data`, a logical vector of
#'   length `nrow(data)`, or a named list with components `name` and `subset`.
#' @param seed Integer seed for deterministic rule evaluation.
#' @param min_n_S Optional exploratory screen for subset size.
#' @param min_varZ_S Optional exploratory screen for within-subset residualized
#'   instrument variance.
#' @param return_all If `TRUE` (default), returns all candidates and diagnostics.
#'   If `FALSE`, returns the first candidate in the user-supplied order. No
#'   automatic subset search or confirmatory selection is performed.
#'
#' @return If `return_all = TRUE`, a data frame with one row per rule and a list
#'   column `subset`. If `return_all = FALSE`, a list with the first supplied
#'   `subset`, `rule`, and `diagnostics`.
#'
#' @details
#' Warning: this function is exploratory. It should not be used for
#' confirmatory BPE inference. Confirmatory BPE requires a pre-specified
#' `bpe_design()` object with a substantive rationale.
#'
#' @examples
#' d <- data.frame(y = rnorm(40), x = rnorm(40), z = rnorm(40),
#'   inactive = rep(c(TRUE, FALSE), each = 20))
#' suppressWarnings(bpe_explore_subsets(
#'   d, y ~ x | z, rules = list(inactive = ~ inactive)))
#'
#' @export
bpe_explore_subsets <- function(data,
                                spec,
                                rules = NULL,
                                seed = 1,
                                min_n_S = NULL,
                                min_varZ_S = NULL,
                                return_all = TRUE) {
  warning(.bpe_exploration_warning_text(), call. = FALSE)
  .bpe_explore_subsets_impl(
    data = data,
    spec = spec,
    rules = rules,
    seed = seed,
    min_n_S = min_n_S,
    min_varZ_S = min_varZ_S,
    return_all = return_all
  )
}

.bpe_explore_subsets_impl <- function(data,
                                      spec,
                                      rules = NULL,
                                      seed = 1,
                                      min_n_S = NULL,
                                      min_varZ_S = NULL,
                                      return_all = TRUE) {
  set.seed(as.integer(seed))

  if (inherits(spec, "formula")) {
    spec <- list(formula = spec)
  }
  if (!is.list(spec) || is.null(spec$formula)) {
    stop("`spec` must be an IV formula or a list containing `formula`.")
  }

  fe <- spec$fe %||% NULL
  fe_engine <- spec$fe_engine %||% "fixest"
  z_names <- spec$z_names %||% NULL
  candidates <- .bpe_normalize_exploratory_rules(rules)

  parsed <- .iv_parse(spec$formula, data, extra_vars = all.vars(fe))

  out_rows <- vector("list", length(candidates))

  for (i in seq_along(candidates)) {
    candidate_i <- candidates[[i]]
    explored_i <- tryCatch(
      .bpe_exploratory_subset(candidate_i, data = data),
      error = function(e) e
    )

    if (inherits(explored_i, "error")) {
      out_rows[[i]] <- data.frame(
        rule = candidate_i$name,
        candidate_type = NA_character_,
        available = FALSE,
        n_S = NA_real_,
        share_S = NA_real_,
        varZ_S = NA_real_,
        F_S = NA_real_,
        screen_n_ok = NA,
        screen_varZ_ok = NA,
        message = conditionMessage(explored_i),
        stringsAsFactors = FALSE
      )
      out_rows[[i]]$subset <- I(list(rep(FALSE, nrow(data))))
      next
    }

    diag_i <- .bpe_subset_diagnostics_from_parsed(
      parsed = parsed,
      data = data,
      subset_idx_full = explored_i$subset,
      fe = fe,
      fe_engine = fe_engine,
      z_names = z_names
    )

    var_min <- if (all(!is.finite(diag_i$varZ_S))) NA_real_ else min(diag_i$varZ_S, na.rm = TRUE)
    msg_i <- if (isTRUE(diag_i$ok)) "" else diag_i$message %||% "Diagnostics unavailable."
    if (length(explored_i$warnings)) {
      msg_i <- paste(c(msg_i, explored_i$warnings), collapse = " | ")
    }

    out_rows[[i]] <- data.frame(
      rule = candidate_i$name,
      candidate_type = explored_i$design$subset_kind,
      available = isTRUE(diag_i$ok),
      n_S = diag_i$n_S,
      share_S = if (is.finite(diag_i$n_S)) diag_i$n_S / length(parsed$y) else NA_real_,
      varZ_S = var_min,
      F_S = as.numeric(diag_i$F_S),
      screen_n_ok = if (is.null(min_n_S)) NA else isTRUE(is.finite(diag_i$n_S) && diag_i$n_S >= min_n_S),
      screen_varZ_ok = if (is.null(min_varZ_S)) NA else isTRUE(is.finite(var_min) && var_min >= min_varZ_S),
      message = msg_i,
      stringsAsFactors = FALSE
    )
    out_rows[[i]]$subset <- I(list(explored_i$subset))
  }

  out <- do.call(rbind, out_rows)

  if (isTRUE(return_all)) {
    return(out)
  }

  sel <- 1L
  list(
    subset = out$subset[[sel]],
    rule = out$rule[sel],
    diagnostics = out,
    message = out$message[sel] %||% ""
  )
}

.embed_prior_by_names <- function(inst_names, z_names, mu_hat, omega_hat) {
  m <- length(inst_names)
  mu <- rep(0, m)
  Omega <- matrix(0, m, m)
  idx <- match(z_names, inst_names)
  if (any(is.na(idx))) {
    stop("Some z_names were not found in instrument names.")
  }
  mu[idx] <- as.numeric(mu_hat)
  Omega[idx, idx] <- as.matrix(omega_hat)
  names(mu) <- inst_names
  dimnames(Omega) <- list(inst_names, inst_names)
  list(mu = mu, Omega = Omega, inst_names = inst_names)
}
