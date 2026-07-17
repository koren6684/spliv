# spliv

[![R-CMD-check](https://github.com/koren6684/spliv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/koren6684/spliv/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

`spliv` provides sensitivity analysis for instrumental-variables (IV) designs when exclusion may fail in structured ways, such as in geolocated clustering or very large data panels. It keeps uniform Conley-style uncertainty as a baseline, adds researcher-specified direct-effect patterns for spatial or panel data, traces uncertainty over sensitivity paths, and supports confirmatory Beyond Plausibly Exogenous (BPE) designs based on a pre-specified instrument-inactive subset. The package does not discover an unrestricted direct-effect field: the pattern or BPE design must be justified by the researcher.

## Installation

```r
remotes::install_github("koren6684/spliv")
```

## A small synthetic example

The examples below are self-contained and use no empirical data.

```r
library(spliv)

set.seed(42)
n <- 240
z <- rnorm(n)
w <- rnorm(n)
inactive <- seq_len(n) <= n / 2
x <- ifelse(inactive, 0, 1) * z + 0.4 * w + rnorm(n)
y <- 1.2 * x + 0.25 * w + 0.15 * z + rnorm(n)
d <- data.frame(y, x, z, w, inactive)
f <- y ~ x + w | z + w

baseline <- spliv(f, d, method = "uci", delta = 0, vcov = "hc1")
baseline$estimates
```

## Uniform UCI sensitivity

Uniform UCI varies the excluded instrument's direct effect over a bounded interval.

```r
uniform <- spliv(
  f, d, method = "uci", delta = 0.20, vcov = "hc1",
  grid = list(steps = 11)
)
uniform$estimates
```

## Patterned UCI and LTZ sensitivity

Use a theory-motivated `spliv_pattern()` to allow direct effects to vary with an observed exposure.

```r
pattern <- spliv_pattern(
  name = "Exposure pattern",
  pattern = ~ w,
  rationale = "The alternative channel is expected to be stronger at higher exposure.",
  variables_used = "w",
  pattern_type = "theory_defined",
  normalize = "max_abs"
)

patterned_uci <- spliv(
  f, d, method = "uci", delta = 0.20, vcov = "hc1",
  violation_pattern = pattern, grid = list(steps = 11)
)
patterned_ltz <- spliv(
  f, d, method = "ltz", delta = 0.20, vcov = "hc1",
  violation_pattern = pattern
)
```

## Sensitivity paths and tipping points

```r
path <- spliv_sensitivity_path(
  f, d, method = "uci", delta_grid = seq(0, 0.30, by = 0.05),
  vcov = "hc1", violation_pattern = pattern
)
head(path)
spliv_tipping_point(path)
plot(path, term = "x")
```

## Confirmatory BPE

BPE begins with an outcome-independent, pre-specified design and validates it before estimation.

```r
design <- bpe_design(
  name = "Theory-defined inactive subset",
  subset = ~ inactive,
  rationale = "The treatment channel is absent in the inactive subset.",
  variables_used = "inactive",
  subset_type = "theory_defined",
  pre_specified = TRUE,
  transportability_rationale = "The subset direct effect is informative for the target sample."
)

validation <- bpe_validate_design(
  f, d, design = design, vcov = "hc1",
  bpe_min_n_S = 40,
  bpe_equiv_margin = 0.25 * sd(resid(lm(x ~ w)))
)
validation[c("n_S", "equivalence_passed", "eligibility_passed")]

# This is a scale-aware illustrative margin for the synthetic example. In a
# substantive analysis, pre-specify the margin from the scientific design; do
# not tune it to make BPE pass.
bpe_margin <- 0.25 * sd(resid(lm(x ~ w)))
bpe_fit <- spliv(
  f, d, method = "bpe", bpe_design = design,
  vcov = "hc1", bpe_min_n_S = 40, bpe_equiv_margin = bpe_margin
)
bpe_fit$estimates
```

`bpe_explore_subsets()` and the deprecated `bpe_find_subset()` are exploratory diagnostics. Searching across subgroups and reporting the first passing rule is **not** confirmatory BPE. Confirmatory BPE requires a pre-specified `bpe_design()` with a substantive rationale and transportability statement.

## Further resources

- Reproducibility repository: [koren6684/spliv-reproducibility](https://github.com/koren6684/spliv-reproducibility) (prospective URL)
- Citation: see [`CITATION.cff`](CITATION.cff) and `citation("spliv")`

The public reproducibility repository contains synthetic simulation workflows and data-gated Koren (2018) and Lelkes, Sood, and Iyengar (2017) replication runners. Restricted third-party data are not redistributed.
