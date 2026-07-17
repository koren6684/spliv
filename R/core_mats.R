.safe_solve <- function(A, B = NULL) {
  if (is.null(B)) {
    return(qr.solve(A, diag(nrow(A))))
  }
  qr.solve(A, B)
}

.drop_intercept_and_constants <- function(M, tol = 1e-10) {
  if (is.null(M)) {
    return(NULL)
  }
  M <- as.matrix(M)
  keep <- rep(TRUE, ncol(M))
  if (!is.null(colnames(M))) {
    keep[colnames(M) == "(Intercept)"] <- FALSE
  }
  vars <- apply(M, 2, function(x) stats::var(as.numeric(x)))
  keep[is.na(vars) | vars < tol] <- FALSE
  if (!any(keep)) {
    return(matrix(0, nrow = nrow(M), ncol = 0))
  }
  M[, keep, drop = FALSE]
}

.iv_2sls_mats <- function(y, X, Z,
                          vcov = c("hc1", "hc0", "iid", "cluster"),
                          cluster_id = NULL) {
  vcov <- match.arg(tolower(vcov), c("hc1", "hc0", "iid", "cluster"))

  y <- as.numeric(y)
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  n <- length(y)
  k <- ncol(X)

  ZZ <- crossprod(Z)
  ZZinv <- .safe_solve(ZZ)
  XZ <- crossprod(X, Z)
  ZY <- crossprod(Z, y)
  XPZX <- XZ %*% ZZinv %*% t(XZ)

  beta <- as.numeric(.safe_solve(XPZX, XZ %*% ZZinv %*% ZY))
  u <- as.numeric(y - X %*% beta)
  bread <- .safe_solve(XPZX)

  if (vcov == "iid") {
    s2 <- sum(u^2) / max(1, n - k)
    V <- s2 * bread
    return(list(beta = beta, vcov = V, resid = u, XZ = XZ, ZZinv = ZZinv, XPZX = XPZX))
  }

  if (vcov %in% c("hc0", "hc1")) {
    Zu <- Z * u
    meat_Z <- crossprod(Zu)
    if (vcov == "hc1") {
      meat_Z <- (n / max(1, n - k)) * meat_Z
    }
    middle <- XZ %*% ZZinv %*% meat_Z %*% ZZinv %*% t(XZ)
    V <- bread %*% middle %*% bread
    return(list(beta = beta, vcov = V, resid = u, XZ = XZ, ZZinv = ZZinv, XPZX = XPZX))
  }

  if (is.null(cluster_id)) {
    stop("vcov='cluster' requires cluster_id.")
  }
  cluster_id <- as.factor(cluster_id)
  G <- nlevels(cluster_id)
  if (G < 2) {
    stop("Need at least 2 clusters for cluster vcov.")
  }

  meat_Z <- matrix(0, ncol(Z), ncol(Z))
  for (g in levels(cluster_id)) {
    idx <- which(cluster_id == g)
    Zg <- Z[idx, , drop = FALSE]
    ug <- u[idx]
    zg_ug <- crossprod(Zg, ug)
    meat_Z <- meat_Z + zg_ug %*% t(zg_ug)
  }

  meat_Z <- (G / (G - 1)) * ((n - 1) / max(1, n - k)) * meat_Z
  middle <- XZ %*% ZZinv %*% meat_Z %*% ZZinv %*% t(XZ)
  V <- bread %*% middle %*% bread
  list(beta = beta, vcov = V, resid = u, XZ = XZ, ZZinv = ZZinv, XPZX = XPZX)
}

sp_ltz_mats <- function(y, X, Z, mu, Omega,
                        level = 0.95,
                        vcov = c("hc1", "hc0", "iid", "cluster"),
                        cluster_id = NULL,
                        coef_names = colnames(X),
                        inst_names = colnames(Z),
                        direct_effect = NULL,
                        direct_effect_names = NULL) {
  fit <- .iv_2sls_mats(y, X, Z, vcov = vcov, cluster_id = cluster_id)

  G <- if (is.null(direct_effect)) {
    as.matrix(Z)
  } else {
    as.matrix(direct_effect)
  }
  if (nrow(G) != nrow(Z)) {
    stop("`direct_effect` must have the same number of rows as `Z`.")
  }

  if (is.null(direct_effect_names)) {
    direct_effect_names <- colnames(G)
  }
  if (is.null(direct_effect_names)) {
    direct_effect_names <- paste0("g", seq_len(ncol(G)))
  }

  mu <- as.numeric(mu)
  if (length(mu) != ncol(G)) {
    stop("`mu` length must match the number of direct-effect columns.")
  }
  if (!is.matrix(Omega)) {
    stop("`Omega` must be a matrix.")
  }
  if (nrow(Omega) != ncol(Omega) || nrow(Omega) != length(mu)) {
    stop("`Omega` must be square with dimensions length(mu) x length(mu).")
  }
  Omega <- as.matrix(Omega)

  ZG <- crossprod(Z, G)
  A <- .safe_solve(fit$XPZX, fit$XZ %*% fit$ZZinv %*% ZG)
  b_adj <- fit$beta - as.numeric(A %*% mu)
  A_mu <- as.numeric(A %*% mu)
  beta_shift <- as.numeric(fit$beta - b_adj)
  V_adj <- fit$vcov + A %*% Omega %*% t(A)

  z <- stats::qnorm(1 - (1 - level) / 2)
  se <- sqrt(pmax(0, diag(V_adj)))
  out <- data.frame(
    term = coef_names,
    estimate = b_adj,
    std.error = se,
    conf.low = b_adj - z * se,
    conf.high = b_adj + z * se,
    row.names = NULL
  )

  list(
    table = out,
    beta = b_adj,
    vcov = V_adj,
    beta_iv = fit$beta,
    iv_beta = fit$beta,
    iv_vcov = fit$vcov,
    mu_used = mu,
    A = A,
    A_mu = A_mu,
    beta_shift = beta_shift,
    inst_names = inst_names,
    direct_effect_names = direct_effect_names
  )
}

conley_ltz_mats <- sp_ltz_mats

.make_gamma_grid <- function(gmin, gmax, steps) {
  p <- length(gmin)
  grids <- lapply(seq_len(p), function(j) seq(gmin[j], gmax[j], length.out = steps[j]))
  as.matrix(expand.grid(grids, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
}

sp_uci_mats <- function(y, X, Z, inst_idx, gamma_grid,
                        level = 0.95,
                        vcov = c("hc1", "hc0", "iid", "cluster"),
                        cluster_id = NULL,
                        coef_names = colnames(X),
                        direct_effect = NULL,
                        direct_effect_names = NULL) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  k <- ncol(X)
  lo <- rep(Inf, k)
  hi <- rep(-Inf, k)
  zcrit <- stats::qnorm(1 - (1 - level) / 2)

  G <- if (is.null(direct_effect)) {
    inst_idx <- as.integer(inst_idx)
    Z[, inst_idx, drop = FALSE]
  } else {
    as.matrix(direct_effect)
  }
  if (nrow(G) != nrow(Z)) {
    stop("`direct_effect` must have the same number of rows as `Z`.")
  }
  if (ncol(gamma_grid) != ncol(G)) {
    stop("`gamma_grid` must have one column per direct-effect regressor.")
  }
  if (is.null(direct_effect_names)) {
    direct_effect_names <- colnames(G)
  }

  for (r in seq_len(nrow(gamma_grid))) {
    gamma <- as.numeric(gamma_grid[r, ])
    y_adj <- y - as.numeric(G %*% gamma)
    fit <- .iv_2sls_mats(y_adj, X, Z, vcov = vcov, cluster_id = cluster_id)
    se <- sqrt(pmax(0, diag(fit$vcov)))
    lo <- pmin(lo, fit$beta - zcrit * se)
    hi <- pmax(hi, fit$beta + zcrit * se)
  }

  data.frame(
    term = coef_names,
    conf.low = lo,
    conf.high = hi,
    row.names = NULL
  )
}

conley_uci_mats <- sp_uci_mats

#' Local-to-Zero Inference for Plausibly Exogenous IV
#'
#' Formula-level wrapper over the matrix core LTZ implementation.
#'
#' @param formula IV formula `y ~ X | Z`.
#' @param data Data frame.
#' @param omega Prior covariance matrix over instrument direct effects.
#' @param mu Prior mean vector over instrument direct effects.
#' @param level Confidence level.
#' @param vcov One of `"iid"`, `"hc0"`, `"hc1"`, `"cluster"`.
#' @param cluster Cluster ids when `vcov = "cluster"`.
#'
#' @return Data frame with adjusted estimates and confidence intervals.
#' @examples
#' set.seed(3)
#' d <- data.frame(y = rnorm(80), x = rnorm(80), z = rnorm(80))
#' prior <- conley_prior_ltz(y ~ x | z, d, inst_vary = "z", sd = 0.1)
#' sp_ltz(y ~ x | z, d, prior$omega, prior$mu)
#' @export
sp_ltz <- function(formula, data, omega, mu,
                       level = 0.95,
                       vcov = c("hc1", "hc0", "iid", "cluster"),
                       cluster = NULL) {
  vcov <- match.arg(tolower(vcov), c("hc1", "hc0", "iid", "cluster"))
  parsed <- .iv_parse(formula, data)
  cluster_id <- NULL
  if (vcov == "cluster") {
    if (is.null(cluster)) {
      stop("vcov='cluster' requires cluster.")
    }
    cluster_id <- as.vector(cluster)[parsed$keep]
  }

  out <- sp_ltz_mats(
    y = parsed$y,
    X = parsed$X,
    Z = parsed$Z,
    mu = mu,
    Omega = omega,
    level = level,
    vcov = vcov,
    cluster_id = cluster_id,
    coef_names = parsed$coef_names,
    inst_names = parsed$inst_names
  )
  out$table
}

#' @rdname sp_ltz
#' @export
conley_ltz <- function(formula, data, omega, mu,
                       level = 0.95,
                       vcov = c("hc1", "hc0", "iid", "cluster"),
                       cluster = NULL) {
  sp_ltz(
    formula = formula,
    data = data,
    omega = omega,
    mu = mu,
    level = level,
    vcov = vcov,
    cluster = cluster
  )
}

#' Union of Confidence Intervals for Plausibly Exogenous IV
#'
#' Formula-level wrapper over the matrix core UCI implementation.
#'
#' @param formula IV formula `y ~ X | Z`.
#' @param data Data frame.
#' @param inst Instrument names to vary.
#' @param gmin Lower bound(s) on gamma.
#' @param gmax Upper bound(s) on gamma.
#' @param grid Grid size per varied instrument.
#' @param level Confidence level.
#' @param vcov One of `"iid"`, `"hc0"`, `"hc1"`, `"cluster"`.
#' @param cluster Cluster ids when `vcov = "cluster"`.
#'
#' @return Data frame with union confidence intervals.
#' @examples
#' set.seed(4)
#' d <- data.frame(y = rnorm(80), x = rnorm(80), z = rnorm(80))
#' sp_uci(y ~ x | z, d, inst = "z", gmin = -0.1, gmax = 0.1, grid = 5)
#' @export
sp_uci <- function(formula, data,
                       inst,
                       gmin,
                       gmax,
                       grid = 21,
                       level = 0.95,
                       vcov = c("hc1", "hc0", "iid", "cluster"),
                       cluster = NULL) {
  vcov <- match.arg(tolower(vcov), c("hc1", "hc0", "iid", "cluster"))
  parsed <- .iv_parse(formula, data)
  if (missing(inst) || length(inst) < 1) {
    stop("Provide inst = c('z1','z2',...) instrument names to vary.")
  }
  if (!all(inst %in% parsed$inst_names)) {
    stop("Some inst names were not found in Z.")
  }
  inst_idx <- match(inst, parsed$inst_names)

  if (length(gmin) == 1) gmin <- rep(gmin, length(inst_idx))
  if (length(gmax) == 1) gmax <- rep(gmax, length(inst_idx))
  if (length(grid) == 1) grid <- rep(grid, length(inst_idx))

  gamma_grid <- .make_gamma_grid(gmin = gmin, gmax = gmax, steps = grid)
  cluster_id <- NULL
  if (vcov == "cluster") {
    if (is.null(cluster)) {
      stop("vcov='cluster' requires cluster.")
    }
    cluster_id <- as.vector(cluster)[parsed$keep]
  }

  sp_uci_mats(
    y = parsed$y,
    X = parsed$X,
    Z = parsed$Z,
    inst_idx = inst_idx,
    gamma_grid = gamma_grid,
    level = level,
    vcov = vcov,
    cluster_id = cluster_id,
    coef_names = parsed$coef_names
  )
}

#' @rdname sp_uci
#' @export
conley_uci <- function(formula, data,
                       inst,
                       gmin,
                       gmax,
                       grid = 21,
                       level = 0.95,
                       vcov = c("hc1", "hc0", "iid", "cluster"),
                       cluster = NULL) {
  sp_uci(
    formula = formula,
    data = data,
    inst = inst,
    gmin = gmin,
    gmax = gmax,
    grid = grid,
    level = level,
    vcov = vcov,
    cluster = cluster
  )
}
