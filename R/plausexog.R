#' @name spliv
#' @section Deprecated Wrappers:
#' `plausexog()` and `plausexog_iv()` are deprecated compatibility aliases for `spliv()`.
#' @aliases plausexog
#' @export
plausexog <- function(formula,
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
  .Deprecated("spliv", package = "spliv")
  spliv(
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

#' @name spliv
#' @aliases plausexog_iv
#' @export
plausexog_iv <- function(...) {
  .Deprecated("spliv", package = "spliv")
  spliv(...)
}
