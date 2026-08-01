# spliv

[![R-CMD-check](https://github.com/koren6684/spliv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/koren6684/spliv/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

`spliv` provides sensitivity analysis for instrumental-variables (IV) designs
when exclusion may fail in structured ways. It supports uniform uncertainty
intervals, researcher-specified direct-effect patterns, sensitivity paths and
tipping points, and confirmatory Beyond Plausibly Exogenous (BPE) designs
based on a pre-specified instrument-inactive subset. The package does not
search for an unrestricted direct-effect field: patterns and BPE designs must
be justified by the researcher.

## Installation

Install the released package from CRAN:

```r
install.packages("spliv")
```

Install the development version from GitHub when you need unreleased changes:

```r
remotes::install_github("koren6684/spliv")
```

## A small synthetic example

The examples below are self-contained and use no empirical or restricted data.

```r
library(spliv)

set.seed(42)
n <- 240
z <- rnorm(n)
w <- rnorm(n)
exposure <- pnorm(w)
inactive <- seq_len(n) <= n / 2
x <- ifelse(inactive, 0, 1) * z + 0.4 * w + rnorm(n)
y <- 1.2 * x + 0.25 * w + 0.15 * exposure * z + rnorm(n)

d <- data.frame(y, x, z, w, exposure, inactive)
f <- y ~ x + w | z + w
```

### Baseline IV estimate

With no method or bound supplied, `spliv()` uses UCI at `delta = 0` and
reproduces the conventional IV confidence interval under strict exclusion.

```r
baseline <- spliv(
  f,
  d,
  vcov = "hc1"
)

baseline$estimates
```

### Uniform UCI sensitivity

Union-of-confidence-intervals (UCI) sensitivity allows the excluded
instrument's direct effect to vary over a bounded interval. On the default
scale, `delta` is an outcome-unit direct effect for a one-residual-SD shift in
the instrument.

```r
uniform <- spliv(
  f,
  d,
  method = "uci",
  delta = 0.20,
  vcov = "hc1",
  grid = list(steps = 11)
)

uniform$estimates
```

### Patterned UCI and LTZ sensitivity

Use a theory-motivated `spliv_pattern()` to allow the possible direct effect to vary with an observed exposure. The package supports both bounded UCI and local-to-zero (LTZ) sensitivity.

```r
pattern <- spliv_pattern(
  name = "Exposure pattern",
  pattern = ~ exposure,
  rationale = "The alternative channel is expected to be stronger at higher exposure.",
  variables_used = "exposure",
  pattern_type = "theory_defined",
  normalize = "max_abs"
)

patterned_uci <- spliv(
  f,
  d,
  method = "uci",
  delta = 0.20,
  vcov = "hc1",
  violation_pattern = pattern,
  grid = list(steps = 11)
)

patterned_ltz <- spliv(
  f,
  d,
  method = "ltz",
  delta = 0.20,
  vcov = "hc1",
  violation_pattern = pattern
)

patterned_uci$estimates
patterned_ltz$estimates
```

### Sensitivity paths and tipping points

Sensitivity paths report how the estimated interval changes over a pre-specified range of direct-effect magnitudes. The tipping point is the first value on the supplied grid at which the interval includes zero.

```r
path <- spliv_sensitivity_path(
  f,
  d,
  method = "uci",
  delta_grid = seq(0, 0.30, by = 0.05),
  vcov = "hc1",
  violation_pattern = pattern
)

head(path)
spliv_tipping_point(path)
plot(path, term = "x")
```

### Confirmatory BPE

Confirmatory Beyond Plausibly Exogenous (BPE) analysis begins with a
pre-specified, outcome-independent instrument-inactive subset and an explicit
transportability rationale. The package validates the proposed design and
estimates the BPE model only when the eligibility diagnostics pass. The default
`sampling` transport carries the estimated reduced-form sampling covariance;
`conservative` adds a pre-specified covariance inflation.

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

# This illustrative margin allows a first-stage effect of 0.25 residual
# treatment SD per one-residual-SD instrument shift. In substantive work,
# pre-specify the margin; do not tune it to make BPE pass.
bpe_margin <- 0.25

validation <- bpe_validate_design(
  f,
  d,
  design = design,
  vcov = "hc1",
  bpe_min_n_S = 40,
  bpe_equiv_margin = bpe_margin
)

validation[c(
  "n_S",
  "equivalence_passed",
  "eligibility_passed"
)]

bpe_fit <- spliv(
  f,
  d,
  method = "bpe",
  bpe_design = design,
  vcov = "hc1",
  bpe_min_n_S = 40,
  bpe_equiv_margin = bpe_margin
)

bpe_fit$estimates
```

`bpe_explore_subsets()` is an exploratory diagnostic. Searching across
subgroups and reporting the first passing rule is **not** confirmatory BPE;
confirmatory BPE requires a pre-specified `bpe_design()` with a substantive
rationale and transportability statement.

## Advanced low-level interface

Advanced users needing lower-level controls can consult the online reference
documentation. Ordinary analyses should normally use the canonical workflow
shown above.

## Further resources

- Package website: <https://koren6684.github.io/spliv/>
- Reproducibility repository: <https://github.com/koren6684/spliv-reproducibility>
- Citation information: run `citation("spliv")` after installation.
