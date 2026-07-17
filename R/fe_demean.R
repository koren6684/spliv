.build_fe_frame <- function(fe_fml, data) {
  if (is.null(fe_fml)) return(NULL)
  if (!inherits(fe_fml, "formula")) {
    stop("`fe` must be a one-sided formula, e.g. ~ gid + ym")
  }
  fe_df <- stats::model.frame(fe_fml, data = data, na.action = stats::na.pass)
  if (ncol(fe_df) < 1) {
    stop("`fe` did not produce any fixed-effect dimensions.")
  }
  fe_df
}

#' Demean with fixest
#'
#' Demeans vectors/matrices by high-dimensional fixed effects using `fixest::demean`.
#'
#' @param y Outcome vector.
#' @param X Regressor matrix.
#' @param Z Instrument matrix.
#' @param W Optional controls matrix.
#' @param fe_fml One-sided FE formula.
#' @param data Data frame aligned with the rows of `y`, `X`, `Z`, `W`.
#'
#' @return List with demeaned `y`, `X`, `Z`, `W`, and FE frame.
#' @examples
#' d <- data.frame(g = factor(rep(1:2, each = 3)))
#' y <- 1:6; X <- matrix(rnorm(6), 6, 1); Z <- matrix(rnorm(6), 6, 1)
#' demean_fixest(y, X, Z, fe_fml = ~ g, data = d)
#' @export
demean_fixest <- function(y, X, Z, W = NULL, fe_fml, data) {
  fe_df <- .build_fe_frame(fe_fml, data)
  if (!requireNamespace("fixest", quietly = TRUE)) {
    stop("Package 'fixest' is required for fe_engine='fixest'.")
  }

  yd <- fixest::demean(as.numeric(y), fe_df)
  Xd <- if (is.null(X) || ncol(X) == 0) X else fixest::demean(as.matrix(X), fe_df)
  Zd <- if (is.null(Z) || ncol(Z) == 0) Z else fixest::demean(as.matrix(Z), fe_df)
  Wd <- if (is.null(W) || ncol(W) == 0) W else fixest::demean(as.matrix(W), fe_df)

  list(y = as.numeric(yd), X = Xd, Z = Zd, W = Wd, fe_df = fe_df)
}

#' Demean with lfe
#'
#' Demeans vectors/matrices by high-dimensional fixed effects using `lfe::demeanlist`.
#'
#' @param y Outcome vector.
#' @param X Regressor matrix.
#' @param Z Instrument matrix.
#' @param W Optional controls matrix.
#' @param fe_list List/data.frame of FE ids.
#'
#' @return List with demeaned `y`, `X`, `Z`, `W`.
#' @examples
#' d <- data.frame(g = factor(rep(1:2, each = 3)))
#' y <- 1:6; X <- matrix(rnorm(6), 6, 1); Z <- matrix(rnorm(6), 6, 1)
#' demean_lfe(y, X, Z, fe_list = d["g"])
#' @export
demean_lfe <- function(y, X, Z, W = NULL, fe_list) {
  if (!requireNamespace("lfe", quietly = TRUE)) {
    stop("Package 'lfe' is required for fe_engine='lfe'.")
  }
  if (is.null(fe_list) || length(fe_list) == 0) {
    stop("`fe_list` must contain at least one fixed-effect id.")
  }
  fl <- as.data.frame(fe_list, stringsAsFactors = FALSE)

  items <- list(y = as.numeric(y), X = as.matrix(X), Z = as.matrix(Z))
  if (!is.null(W)) {
    items$W <- as.matrix(W)
  }
  out <- lfe::demeanlist(items, fl = fl)

  list(
    y = as.numeric(out$y),
    X = out$X,
    Z = out$Z,
    W = out$W %||% W,
    fe_df = fl
  )
}
