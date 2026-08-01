# Sensitivity-path benchmark

`benchmark_sensitivity_paths.R` compares equivalent UCI computations:

1. `spliv_sensitivity_path()`; and
2. a transparent loop that adjusts the outcome at the identical theta grid and
   refits the identical residualized IV design with `fixest::feols()`.

Both implementations use the same normalized pattern, residual-SD-scaled
instrument, theta and delta grids, HC1 covariance convention, and union rule.
The script stops if interval endpoints differ by more than `1e-8`; an unequal
comparison is never timed as a valid speed comparison.

The full design covers sample sizes 1,000, 10,000, and 50,000; multiple
fixed-effect group counts; delta-grid lengths 5, 21, and 41; uniform and
patterned sensitivity; and no-FE and two-way-FE designs. Each method receives a
warm-up run followed by at least five timed iterations. Median time, IQR,
iterations per second, `Rprofmem()` allocated bytes, errors, model dimensions,
software versions, operating system, and seed are recorded.

Run from the package root:

```sh
Rscript dev/benchmarks/benchmark_sensitivity_paths.R
```

The completed full matrix is recorded in
`output/benchmark_results_full.csv`. A representative correctness and timing
smoke run uses the same five-iteration rule:

```sh
SPLIV_BENCHMARK_FAST=TRUE Rscript dev/benchmarks/benchmark_sensitivity_paths.R
```

`SPLIV_BENCHMARK_THETA_STEPS` controls the inner theta grid (default 5) and
`SPLIV_BENCHMARK_ITERATIONS` may increase, but not reduce below five, the
number of timed iterations. Outputs are written under `dev/benchmarks/output/`
and are excluded from the CRAN source package.

No experimental cached-design implementation is included. Adding one safely
would require a separate public contract for reusable parsed/residualized
designs; this release does not expose such an interface.
