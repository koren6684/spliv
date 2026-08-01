approved <- sort(c(
  "spliv", "spliv_pattern", "spliv_eval_pattern",
  "spliv_sensitivity_path", "spliv_tipping_point",
  "bpe_design", "bpe_eval_subset", "bpe_validate_design",
  "bpe_explore_subsets", "sp_ltz", "sp_uci", "sp_prior_ltz",
  "sp_sensitivity_ltz_normal",
  "sp_sensitivity_ltz_uniform01_as_normal",
  "sp_sensitivity_uci_support", "plot_sp_sensitivity"
))
removed <- c(
  "plausexog", "plausexog_iv", "conley_ltz", "conley_uci",
  "conley_prior_ltz", "conley_sensitivity_ltz_normal",
  "conley_sensitivity_uci_support",
  "conley_sensitivity_ltz_uniform01_as_normal",
  "plot_conley_sensitivity", "bpe_find_subset",
  "estimate_gamma_zero_first_stage", "demean_fixest", "demean_lfe",
  "embed_prior_into_full_Z", "iv_inst_names"
)

stopifnot(as.character(packageVersion("spliv")) == "0.2.0")
stopifnot(identical(sort(getNamespaceExports("spliv")), approved))
stopifnot(!any(vapply(removed, exists, logical(1), envir = asNamespace("spliv"), inherits = FALSE)))

set.seed(9001)
n <- 4000
z <- rnorm(n)
w <- rnorm(n)
inactive <- seq_len(n) <= n / 2
x <- ifelse(inactive, 0, 1) * z + 0.4 * w + rnorm(n)
y <- 1.2 * x + 0.25 * w + 0.1 * z + rnorm(n)
d <- data.frame(y, x, z, w, exposure = pnorm(w), inactive)
f <- y ~ x + w | z + w

default_fit <- spliv::spliv(f, d)
explicit_fit <- spliv::spliv(f, d, method = "uci", delta = 0)
stopifnot(identical(class(default_fit), "spliv_fit"))
stopifnot(isTRUE(all.equal(default_fit$estimates, explicit_fit$estimates, tolerance = 1e-12)))

pattern <- spliv::spliv_pattern(
  "Exposure", ~ exposure,
  rationale = "The alternative channel follows exposure."
)
patterned <- spliv::spliv(f, d, method = "uci", delta = 0.1,
                          violation_pattern = pattern, grid = list(steps = 5))
path <- spliv::spliv_sensitivity_path(
  f, d, method = "uci", delta_grid = c(0, 0.05),
  violation_pattern = pattern, grid = list(steps = 5)
)
stopifnot(inherits(patterned, "spliv_fit"), inherits(path, "spliv_sensitivity_path"))

plot_file <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_file)
stopifnot(identical(plot(default_fit), default_fit))
stopifnot(identical(plot(path, term = "x"), path))
grDevices::dev.off()
unlink(plot_file)

valid_design <- spliv::bpe_design(
  "Inactive subset", ~ inactive,
  rationale = "The treatment channel is absent in this pre-specified subset.",
  transportability_rationale = "The direct-effect mechanism is assumed to apply to the target sample."
)
validation <- spliv::bpe_validate_design(
  f, d, valid_design, bpe_min_n_S = 1000, bpe_equiv_margin = 0.25
)
stopifnot(isTRUE(validation$eligibility_passed))
bpe_fit <- spliv::spliv(
  f, d, method = "bpe", bpe_design = valid_design,
  bpe_min_n_S = 1000, bpe_equiv_margin = 0.25
)
stopifnot(inherits(bpe_fit, "spliv_fit"), all(is.finite(bpe_fit$estimates$conf.low)))

invalid_design <- spliv::bpe_design(
  "Missing transport rationale", ~ inactive,
  rationale = "The treatment channel is absent in this pre-specified subset."
)
invalid <- spliv::bpe_validate_design(
  f, d, invalid_design, bpe_min_n_S = 1000, bpe_equiv_margin = 0.25
)
stopifnot(!isTRUE(invalid$eligibility_passed))
withheld <- spliv::spliv(
  f, d, method = "bpe", bpe_design = invalid_design,
  bpe_min_n_S = 1000, bpe_equiv_margin = 0.25,
  bpe_not_applicable = "na"
)
stopifnot(all(is.na(withheld$estimates$estimate)))

cat("clean-room package test: PASS\n")
cat("version:", as.character(packageVersion("spliv")), "\n")
cat("library:", find.package("spliv"), "\n")
cat("exports:", paste(approved, collapse = ", "), "\n")
