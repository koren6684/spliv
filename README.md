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
inactive <- seq_len(n) <= n / 2
x <- ifelse(inactive, 0, 1) * z + 0.4 * w + rnorm(n)
y <- 1.2 * x + 0.25 * w + 0.15 * z + rnorm(n)
d <- data.frame(y, x, z, w, inactive)
f <- y ~ x + w | z + w

fit <- spliv(f, d, method = "ltz", delta = 0, vcov = "hc1")
fit$estimates
plot(fit)
```

## Uniform UCI

Uniform UCI varies an excluded instrument's direct effect over a bounded
interval:

```r
uniform <- spliv(
  f, d, method = "uci", delta = 0.20, vcov = "hc1",
  grid = list(steps = 11)
)
uniform$estimates
```

## Patterned UCI/LTZ

Use a theory-motivated `spliv_pattern()` when the plausible direct effect is
expected to vary with an observed exposure:

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

Trace a pre-specified grid and report the first grid value at which an interval
contains zero:

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

BPE starts with an outcome-independent, pre-specified design and validates it
before estimation. The margin below is a scale-aware, synthetic illustrative
choice; in substantive work, pre-specify it from the scientific design rather
than tuning it after seeing whether BPE passes.

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

bpe_margin <- 0.25 * sd(resid(lm(x ~ w)))
validation <- bpe_validate_design(
  f, d, design = design, vcov = "hc1",
  bpe_min_n_S = 40, bpe_equiv_margin = bpe_margin
)
validation[c("n_S", "equivalence_passed", "eligibility_passed")]

bpe_fit <- spliv(
  f, d, method = "bpe", bpe_design = design,
  vcov = "hc1", bpe_min_n_S = 40, bpe_equiv_margin = bpe_margin
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
