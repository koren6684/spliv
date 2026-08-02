# Sensitivity-path benchmark

`benchmark_sensitivity_paths.R` compares equivalent UCI computations:

1. `spliv_sensitivity_path()`; and
2. a transparent loop that adjusts the outcome at the identical theta grid and
   refits the same residualized IV design with `fixest::feols()`.

Both implementations use the same normalized pattern, residual-SD-scaled
instrument, theta and delta grids, HC1 covariance convention, and UCI union
rule. The benchmark stops if corresponding interval endpoints differ by more
than `1e-8`; a numerically nonequivalent comparison is never treated as a
valid timing comparison.

The optimized SPLIV implementation prepares the invariant components of the
model once per sensitivity path. Formula parsing, complete-case alignment,
pattern evaluation, instrument scaling, fixed-effect residualization, cluster
alignment, and invariant design calculations are then reused across delta
values. The transparent comparison loop instead adjusts the outcome and refits
the model separately at each sensitivity value.

The full benchmark covers sample sizes of 1,000, 10,000, and 50,000; multiple
fixed-effect group-count designs; delta-grid lengths of 5, 21, and 41; uniform
and patterned sensitivity; and specifications with no fixed effects or two-way
fixed effects. Each method receives one warm-up run followed by at least five
timed iterations. The output records median runtime, timing IQR, iterations per
second, `Rprofmem()` allocated bytes, errors, model dimensions, software
versions, operating system, and random seed.

Run the full benchmark from the package root:

```sh
Rscript dev/benchmarks/benchmark_sensitivity_paths.R
