# Benchmark findings

## Status

The complete 72-case equivalent-computation matrix passed on the release
hardware. It covered `n = 1,000`, `10,000`, and `50,000`; two fixed-effect
group-count designs per sample size; 5-, 21-, and 41-value delta grids;
uniform and patterned UCI; and no-FE and two-way-FE specifications. Each
method received one warm-up and five measured iterations. The machine-readable
results are in `output/benchmark_results_full.csv` (SHA-256
`<NEW_SHA256>`).

All 144 method rows completed without error. The maximum absolute disagreement
between corresponding UCI endpoints was `5.24e-14`, well below the
pre-specified `1e-8` numerical-equivalence tolerance.

## Main results

The benchmark was run on Darwin 23.5.0 arm64 with R 4.5.2. The optimization
substantially reduced sensitivity-path runtime while preserving the resulting
estimates and confidence intervals. Relative to the pre-optimization
implementation, the median speedup for `spliv_sensitivity_path()` was `2.84`
times in models without fixed effects and `6.86` times in models with two-way
fixed effects. The repeated adjusted-outcome `fixest` comparator changed by
only about `1.10` times across the same comparisons, indicating that the larger
SPLIV gains reflect the path optimization rather than a general change in
machine performance.

The optimized SPLIV implementation was faster than repeated adjusted-outcome
refitting in every benchmarked design. In the least favorable comparison, a
uniform two-way-fixed-effect design with `n = 50,000`, SPLIV required `0.263`
seconds compared with `0.311` seconds for repeated refitting. Its time ratio
was `0.846`, corresponding to a runtime reduction of approximately `15.4%`.
In the most favorable comparison, a patterned design with `n = 1,000`, SPLIV
required `0.008` seconds compared with `1.253` seconds for repeated refitting.
Its time ratio was `0.0064`, making SPLIV approximately `157` times faster in
that case.

The runtime gains do not imply uniformly lower memory use. In the least
favorable timing comparison, SPLIV allocated approximately `419.6` MB compared
with `336.4` MB for repeated refitting, a memory ratio of `1.247`. In the most
favorable comparison, SPLIV allocated approximately `1.7` MB compared with
`92.8` MB, a ratio of `0.018`. The optimized implementation therefore reduced
runtime throughout the benchmark while memory performance continued to depend
on the particular workload.

## Interpretation and limits

The transparent comparator adjusts the outcome at each sensitivity value and
refits the same residualized IV design using `fixest`. It uses the same
normalized instrument, violation pattern, theta grid, delta grid, HC1
convention, and UCI union rule as SPLIV. The numerical-equivalence checks
confirm that the two approaches produce corresponding confidence-region
endpoints to machine precision.

The optimized SPLIV path now prepares the invariant components of the analysis
once per sensitivity path. Formula parsing, complete-case alignment, pattern
evaluation, instrument scaling, fixed-effect residualization, cluster
alignment, and invariant design calculations are reused across delta values
rather than reconstructed at every step. This change is especially consequential
for two-way-fixed-effect panels, where the median speedup relative to the prior
SPLIV implementation was nearly sevenfold.

These results support two conclusions. First, SPLIV makes patterned
exclusion-restriction sensitivity analysis computationally practical for large
fixed-effect panels. Second, its optimized path was faster than transparent
repeated adjusted-outcome refitting in every workload represented in this
benchmark. The evidence does not establish uniformly lower memory use or
performance superiority for every possible model.

This remains a single-machine benchmark with five measured iterations, a
five-value inner theta grid, one endogenous treatment, one excluded instrument,
one exogenous control, and HC1 inference. Performance may differ on other
hardware, with clustered covariance estimation, additional controls,
multiple model terms, or other workloads outside the current high-level
package contract.
