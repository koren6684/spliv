# Instrument Names for Parsed IV Formula

Returns the instrument-matrix column order used internally.

## Usage

``` r
iv_inst_names(formula, data)
```

## Arguments

- formula:

  IV formula `y ~ X | Z`.

- data:

  Data frame.

## Value

Character vector of instrument names.

## Examples

``` r
d <- data.frame(y = 1:4, x = 1:4, z = 4:1)
iv_inst_names(y ~ x | z, d)
#> [1] "(Intercept)" "z"          
```
