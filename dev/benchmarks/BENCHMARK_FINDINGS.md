# Benchmark findings

## Status

The complete 72-case equivalent-computation matrix passed on the release
hardware. It covered `n = 1,000`, `10,000`, and `50,000`; two fixed-effect
group-count designs per sample size; 5-, 21-, and 41-value delta grids;
uniform and patterned UCI; and no-FE and two-way-FE specifications. Each
method received one warm-up and five measured iterations. The machine-readable
results are in `output/benchmark_results_full.csv` (SHA-256
`536eadbbe298486cf57fe27937338db9a366f423311a0691e604bcbafa274173`).

All 144 method rows completed without error. Maximum absolute disagreement in
a union endpoint was `5.24e-14`, well below the pre-specified `1e-8`
equivalence tolerance.

## Main results

The comparison was run on Darwin 23.5.0 arm64 with R 4.5.2 and spliv 0.2.0.
Across all cases, the median SPLIV/reference elapsed-time ratio was `0.283`
(range `0.069` to `1.259`), and the median allocation ratio was `0.685`
(range `0.339` to `1.173`). Median elapsed times across the case matrix were
`0.222` seconds for SPLIV and `0.864` seconds for transparent repeated
adjusted-outcome `fixest` refits.

| n | SPLIV median (s) | Refit median (s) | Median time ratio | SPLIV median bytes | Refit median bytes | Median byte ratio |
|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 0.0830 | 0.7335 | 0.175 | 25,083,900 | 40,768,248 | 0.660 |
| 10,000 | 0.2065 | 0.8640 | 0.337 | 251,276,748 | 354,508,248 | 0.772 |
| 50,000 | 0.9630 | 1.6115 | 0.773 | 1,240,297,248 | 1,748,908,248 | 0.776 |

The aggregate result hides an important fixed-effect crossover. With no fixed
effects, the median SPLIV/reference time ratios were `0.080`, `0.158`, and
`0.346` at `n = 1,000`, `10,000`, and `50,000`; median allocation ratios were
`0.363`, `0.412`, and `0.411`. With two-way fixed effects, the corresponding
time ratios were `0.267`, `0.516`, and `1.196`, while allocation ratios were
`0.959`, `1.135`, and `1.145`. Thus SPLIV was faster in every no-FE case and
in the smaller two-way-FE cases, but it was 10.7% to 25.9% slower and allocated
11.7% to 17.3% more memory in every `n = 50,000` two-way-FE case.

Grid length scaled close to linearly for both implementations: median SPLIV
times were `0.075`, `0.304`, and `0.595` seconds for 5, 21, and 41 delta
values; comparator medians were `0.215`, `0.864`, and `1.692` seconds.
Patterned cases were modestly more expensive than uniform cases for SPLIV
(`0.251` versus `0.222` seconds at the cross-case median). The largest timing
IQR recorded was `0.097` seconds.

## Interpretation and limits

The transparent comparator refits the same residualized IV design at every
theta using the identical normalized instrument, pattern, theta grid, delta
grid, HC1 convention, and UCI union rule. SPLIV reuses its matrix design within
each UCI union, explaining its substantial advantage without fixed effects.
SPLIV still repeats formula parsing and fixed-effect residualization across
delta values; that repeated work dominates enough at `n = 50,000` for the
two-way-FE implementation to become slower and more allocation-intensive than
the reference.

This is a single-machine, five-iteration benchmark with a five-value inner
theta grid and one treatment, one instrument, and one control. It does not
establish performance on other hardware, richer formulas, or other covariance
estimators. No experimental cached-across-delta implementation was added, so
the potential benefit of such caching remains unknown. The evidence supports
a qualified workload-specific result, not an unconditional claim that SPLIV
is fast or computationally efficient.
