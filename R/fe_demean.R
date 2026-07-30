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

.demean_fixest <- function(y, X, Z, W = NULL, fe_fml, data) {
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

.demean_lfe <- function(y, X, Z, W = NULL, fe_list) {
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
