.trim_ws <- function(x) gsub("^\\s+|\\s+$", "", x)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.iv_parse <- function(formula, data, extra_vars = character(0)) {
  ftxt <- paste(deparse(formula), collapse = " ")
  parts <- strsplit(ftxt, "\\|", fixed = FALSE)[[1]]
  if (length(parts) != 2) {
    stop("Formula must be of the form: y ~ X | Z")
  }

  main_f <- stats::as.formula(.trim_ws(parts[1]))
  inst_f <- stats::as.formula(paste("~", .trim_ws(parts[2])))

  yvar <- all.vars(main_f)[1]
  vars_all <- unique(c(all.vars(main_f), all.vars(inst_f), extra_vars))
  missing_vars <- setdiff(vars_all, names(data))
  if (length(missing_vars) > 0) {
    stop(
      "These variables are referenced in the formula but not found in `data`: ",
      paste(missing_vars, collapse = ", ")
    )
  }

  mf_all <- data[, vars_all, drop = FALSE]
  keep <- stats::complete.cases(mf_all)
  if (!any(keep)) {
    stop("After listwise deletion (complete.cases), no rows remain.")
  }
  mf_all <- mf_all[keep, , drop = FALSE]

  y <- mf_all[[yvar]]
  X <- stats::model.matrix(stats::delete.response(stats::terms(main_f)), data = mf_all)
  Z <- stats::model.matrix(stats::terms(inst_f), data = mf_all)

  if (!"(Intercept)" %in% colnames(Z)) {
    Z <- cbind("(Intercept)" = 1, Z)
  }

  list(
    y = as.numeric(y),
    X = as.matrix(X),
    Z = as.matrix(Z),
    n = length(y),
    k = ncol(X),
    m = ncol(Z),
    coef_names = colnames(X),
    inst_names = colnames(Z),
    keep = keep,
    data_cc = mf_all,
    main_formula = main_f,
    inst_formula = inst_f
  )
}

#' Instrument Names for Parsed IV Formula
#'
#' Returns the instrument-matrix column order used internally.
#'
#' @param formula IV formula `y ~ X | Z`.
#' @param data Data frame.
#'
#' @return Character vector of instrument names.
#' @examples
#' d <- data.frame(y = 1:4, x = 1:4, z = 4:1)
#' iv_inst_names(y ~ x | z, d)
#' @export
iv_inst_names <- function(formula, data) {
  .iv_parse(formula, data)$inst_names
}
