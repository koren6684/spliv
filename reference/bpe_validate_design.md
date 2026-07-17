# Validate a Confirmatory BPE Design

Runs the confirmatory BPE eligibility diagnostics for a pre-specified
design without fitting the final SPLIV model.

## Usage

``` r
bpe_validate_design(
  formula,
  data,
  design,
  fe = NULL,
  fe_engine = c("fixest", "lfe"),
  vcov = c("iid", "hc1", "cluster"),
  cluster = NULL,
  z_names = NULL,
  bpe_min_n_S = 2000,
  bpe_min_clusters_S = 30,
  bpe_min_varZ_S = 1e-06,
  bpe_equiv_margin,
  bpe_equiv_level = 0.95,
  bpe_transport = c("none", "sampling", "conservative"),
  bpe_transport_kappa = 0,
  bpe_kappa = 1,
  scale_instrument = c("residual_sd", "none")
)
```

## Arguments

- formula:

  IV formula `y ~ X | Z`.

- data:

  Data frame.

- design:

  A
  [`bpe_design()`](https://koren6684.github.io/spliv/reference/bpe_design.md)
  object.

- fe:

  Optional one-sided formula of fixed effects.

- fe_engine:

  FE demeaning engine, one of `"fixest"` or `"lfe"`.

- vcov:

  One of `"iid"`, `"hc1"`, or `"cluster"`.

- cluster:

  Cluster ids or one-sided formula when `vcov = "cluster"`.

- z_names:

  Optional instrument name for BPE. Confirmatory BPE currently supports
  exactly one instrument.

- bpe_min_n_S:

  Minimum subset size threshold.

- bpe_min_clusters_S:

  Minimum number of clusters required in the subset when clustered
  covariance is used.

- bpe_min_varZ_S:

  Minimum residualized instrument variance required in the subset.

- bpe_equiv_margin:

  Researcher-specified equivalence margin for the first-stage
  coefficient. Eligibility is based on the first-stage equivalence
  interval, not on the first-stage F-statistic.

- bpe_equiv_level:

  Confidence level used for the first-stage equivalence interval.

- bpe_transport:

  One of `"none"`, `"sampling"`, or `"conservative"`. Transportability
  is an assumption reflected in the reported covariance; it is not
  established by the subset itself.

- bpe_transport_kappa:

  Non-negative scalar controlling the conservative transport covariance
  inflation.

- bpe_kappa:

  Positive scalar multiplier applied to the transported BPE covariance
  before it is embedded into the LTZ prior.

- scale_instrument:

  One of `"residual_sd"` or `"none"`.

## Value

An object of class `"spliv_bpe_validation"` containing design metadata,
subset diagnostics, first-stage equivalence diagnostics, reduced- form
direct-effect estimates, and covariance components for confirmatory BPE.

## Examples

``` r
set.seed(2)
d <- data.frame(
  y = rnorm(80), x = rnorm(80), z = rnorm(80),
  inactive = rep(c(TRUE, FALSE), each = 40)
)
design <- bpe_design("Inactive", ~ inactive,
  rationale = "The treatment channel is absent.")
bpe_validate_design(y ~ x | z, d, design,
  bpe_min_n_S = 20, bpe_equiv_margin = 1)
#> $design_name
#> [1] "Inactive"
#> 
#> $rationale
#> [1] "The treatment channel is absent."
#> 
#> $subset_type
#> NULL
#> 
#> $variables_used
#> [1] "inactive"
#> 
#> $pre_specified
#> [1] TRUE
#> 
#> $transportability_rationale
#> NULL
#> 
#> $notes
#> NULL
#> 
#> $n_S
#> [1] 40
#> 
#> $share_S
#> [1] 0.5
#> 
#> $G_S
#> NULL
#> 
#> $varZ_S
#>         z 
#> 0.8737843 
#> 
#> $residualized_instrument_sd_S
#>         z 
#> 0.9347643 
#> 
#> $residualized_instrument_sd
#>        z 
#> 1.031161 
#> 
#> $residualized_treatment_sd_S
#>        x 
#> 1.150203 
#> 
#> $first_stage_coefficient
#>         z 
#> 0.2973111 
#> 
#> $first_stage_se
#>         z 
#> 0.1936949 
#> 
#> $first_stage_ci
#>         lower     upper
#> z -0.08232386 0.6769461
#> 
#> $first_stage_f_statistic
#>        z 
#> 2.356058 
#> 
#> $first_stage_f_type
#> [1] "conventional_ols_diagnostic"
#> 
#> $first_stage_effect_one_residual_sd_Z
#>         z 
#> 0.2779158 
#> 
#> $standardized_first_stage_effect
#>         z 
#> 0.2416232 
#> 
#> $equivalence_margin
#> z 
#> 1 
#> 
#> $equivalence_level
#> [1] 0.95
#> 
#> $equivalence_passed
#> [1] TRUE
#> 
#> $eligibility_passed
#> [1] TRUE
#> 
#> $eligibility_checks
#> $eligibility_checks$pre_specified
#> [1] TRUE
#> 
#> $eligibility_checks$rationale
#> [1] TRUE
#> 
#> $eligibility_checks$minimum_n
#> [1] TRUE
#> 
#> $eligibility_checks$minimum_clusters
#> [1] TRUE
#> 
#> $eligibility_checks$residual_variation
#> [1] TRUE
#> 
#> $eligibility_checks$equivalence
#> [1] TRUE
#> 
#> 
#> $reduced_form_direct_effect
#>          z 
#> -0.2241848 
#> 
#> $reduced_form_direct_effect_cov
#>            z
#> z 0.03586035
#> 
#> $reduced_form_sampling_cov
#>            z
#> z 0.03586035
#> 
#> $transport_covariance
#>            z
#> z 0.03586035
#> 
#> $transport_mode
#> [1] "none"
#> 
#> $transport_uncertainty_inflation
#> [1] 1
#> 
#> $prior_mu_sub
#>          z 
#> -0.2241848 
#> 
#> $prior_Omega_sub
#>            z
#> z 0.03586035
#> 
#> $prior_mu_full
#> (Intercept)           z 
#>   0.0000000  -0.2241848 
#> 
#> $prior_Omega_full
#>             (Intercept)          z
#> (Intercept)           0 0.00000000
#> z                     0 0.03586035
#> 
#> $subset_idx_full
#>  [1]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#> [13]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#> [25]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#> [37]  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [49] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [61] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [73] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> attr(,"bpe_warnings")
#> character(0)
#> attr(,"bpe_na_share")
#> [1] 0
#> 
#> $design_audit
#> $design_audit$variables_used
#> [1] "inactive"
#> 
#> $design_audit$uses_outcome
#> [1] FALSE
#> 
#> $design_audit$uses_endogenous
#> character(0)
#> 
#> $design_audit$uses_instrument
#> character(0)
#> 
#> $design_audit$diagnostic_hit
#> [1] FALSE
#> 
#> $design_audit$warnings
#> character(0)
#> 
#> $design_audit$expression
#> [1] "~inactive"
#> 
#> $design_audit$instrument
#> [1] "z"
#> 
#> 
#> $warnings
#> character(0)
#> 
#> $message
#> [1] ""
#> 
#> $first_stage_target
#> [1] "x"
#> 
#> $instrument
#> [1] "z"
#> 
#> $scale_instrument
#> [1] "residual_sd"
#> 
#> attr(,"class")
#> [1] "spliv_bpe_validation"
```
