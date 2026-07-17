# Plot Patterned Sensitivity Output or a Fitted Object

Plot Patterned Sensitivity Output or a Fitted Object

## Usage

``` r
plot_sp_sensitivity(
  df_plot,
  ylab = "Effect (beta)",
  main = "Plausibly Exogenous IV Sensitivity",
  ...
)

plot_conley_sensitivity(
  df_plot,
  ylab = "Effect (beta)",
  main = "Plausibly Exogenous IV Sensitivity",
  ...
)
```

## Arguments

- df_plot:

  Data frame returned by a sensitivity helper or `plausexog_fit` object.

- ylab:

  Y-axis label.

- main:

  Plot title.

- ...:

  Unused.

## Value

Invisibly returns the plotted input.

## Examples

``` r
set.seed(7)
d <- data.frame(y = rnorm(80), x = rnorm(80), z = rnorm(80))
p <- spliv_sensitivity_path(y ~ x | z, d, method = "uci", delta_grid = c(0, 0.1))
plot_sp_sensitivity(p, term = "x")
```
