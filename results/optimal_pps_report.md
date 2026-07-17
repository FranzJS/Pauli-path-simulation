# Whole-frontier variance-optimal PPS benchmark

## Algorithm

The deterministic L1-mass BFS pass supplies the retained-support schedule
`K_1, K_2, ...`. At each truncation, the randomized pass:

1. uses the corresponding L1-derived budget `K_i`;
2. solves for `tau` such that `sum_j min(1, abs(c_j) / tau) = K_i`;
3. sets `pi_j = min(1, abs(c_j) / tau)`;
4. uses randomized-order systematic PPS to select exactly `K_i` distinct
   coordinates with marginal inclusion probabilities `pi_j`; and
5. stores a selected coefficient as `c_j / pi_j`.

This is conditionally unbiased and minimizes the one-truncation expected
squared L2 reconstruction error among coordinatewise HT sparsifiers with
expected support `K_i`. Probability-one coordinates are retained
deterministically, so the algorithm derives its heavy-hitter count instead of
requiring an `alpha` parameter.

## Validation

The unit tests check:

- exact support size for every seed;
- deterministic retention of all saturated (`pi_i = 1`) coordinates;
- the zero-mass uniform fallback;
- empirical marginal inclusion probabilities; and
- empirical coordinatewise HT unbiasedness.

The statistical test uses coefficients `(1, 2, 3, 4)` with `K = 2`, for which
the exact probabilities are `(0.2, 0.4, 0.6, 0.8)`. Every selected adjusted
coefficient equals `5`, and the empirical coordinate means reproduce the
inputs. The complete C++ test binary passes.

## Single-pass comparison

These measurements use a native GCC release build. Runtime is the median of
three fresh-process repetitions; RSS is the maximum across those repetitions.

| Case | Method | Absolute error | Runtime (s) | Runtime / BFS | Peak RSS (MB) | RSS / BFS |
|---|---|---:|---:|---:|---:|---:|
| Clifford+T | L1 BFS | 6.905e-4 | 0.658 | 1.00 | 68.59 | 1.00 |
| Clifford+T | optimal PPS HT | 1.763e-3 | 0.819 | 1.24 | 87.24 | 1.27 |
| Ising | L1 BFS | 4.503e-4 | 2.861 | 1.00 | 64.28 | 1.00 |
| Ising | optimal PPS HT | 1.967e-3 | 5.882 | 2.06 | 149.76 | 2.33 |
| Noisy Clifford+T | L1 BFS | 4.616e-4 | 0.00860 | 1.00 | 13.52 | 1.00 |
| Noisy Clifford+T | optimal PPS HT | 4.606e-2 | 0.01488 | 1.73 | 13.49 | 1.00 |

## Repeated-pass error

| Case | Passes | Running-mean absolute error | Per-pass SD | Empirical SE | Mean runtime / pass (s) | Maximum HT multiplier |
|---|---:|---:|---:|---:|---:|---:|
| Clifford+T | 10 | 5.835e-4 | 3.734e-3 | 1.181e-3 | 0.855 | 1.53e3 |
| Clifford+T | 20 | 4.667e-4 | 3.849e-3 | 8.606e-4 | 0.823 | 2.40e3 |
| Clifford+T | 50 | 7.444e-4 | 4.804e-3 | 6.793e-4 | 0.808 | 2.44e3 |
| Clifford+T | 100 | 5.756e-4 | 4.286e-3 | 4.286e-4 | 0.822 | 2.16e4 |
| Ising | 10 | 3.254e-4 | 9.587e-4 | 3.032e-4 | 5.494 | 1.39e4 |
| Ising | 20 | 3.211e-5 | 1.462e-3 | 3.269e-4 | 5.416 | 1.39e4 |
| Ising | 50 | 3.390e-4 | 1.824e-3 | 2.579e-4 | 5.341 | 1.69e4 |
| Ising | 100 | 3.360e-5 | 1.722e-3 | 1.722e-4 | 5.367 | 2.19e4 |
| Noisy Clifford+T | 10 | 4.676e-3 | 3.241e-2 | 1.025e-2 | 0.0156 | 87.3 |
| Noisy Clifford+T | 20 | 1.965e-3 | 2.862e-2 | 6.399e-3 | 0.0158 | 87.3 |
| Noisy Clifford+T | 50 | 2.054e-3 | 2.597e-2 | 3.673e-3 | 0.0155 | 183 |
| Noisy Clifford+T | 100 | 1.986e-3 | 2.468e-2 | 2.468e-3 | 0.0151 | 214 |

## Interpretation

The method behaves as intended mathematically: its support is exactly `K_i`,
its saturated coordinates are deterministic, and the HT estimator is
coordinatewise unbiased. It does not, however, make a pass as cheap as L1 BFS
on these circuits. Rescaled coefficients change later merging and cancellation,
so the pre-truncation frontier can be much larger even though every
post-truncation support matches the baseline schedule. This is most pronounced
for Ising.

At 100 passes, the deviations from the reference are 1.34 empirical standard
errors for Clifford+T, 0.20 for Ising, and 0.80 for noisy Clifford+T. All three
are therefore statistically consistent with unbiased convergence. The absolute
error is not monotone in the number of passes: rare large HT multipliers and
ordinary sampling fluctuations can move the running mean away from the
reference temporarily. In particular, the unusually accurate 20-pass Ising
mean should not be read as a permanently attained error level.

Clifford+T roughly matches deterministic BFS accuracy by 100 passes, and the
Ising mean is substantially more accurate in this realization. Noisy
Clifford+T remains high variance and has a larger realized error than
deterministic BFS, although its 100-pass error is below its empirical standard
error. A noise-aware inclusion score, for example coefficient magnitude times
estimated remaining damping, is the natural next variance-reduction experiment.

The reported execution is sequential to keep runtime and RSS measurements
comparable. The randomized passes share no mutable state and can be distributed
across workers with only a final reduction of estimates and diagnostics.

## Pass-count implication

Using the 100-pass empirical standard deviation `s_100`, a rough estimate of
the number of passes needed to make the Monte Carlo standard error equal to the
observed deterministic BFS absolute error is

`R_match = (s_100 / BFS_error)^2`.

| Case | `s_100` | BFS absolute error | Estimated `R_match` |
|---|---:|---:|---:|
| Clifford+T | 4.286e-3 | 6.905e-4 | 39 |
| Ising | 1.722e-3 | 4.503e-4 | 15 |
| Noisy Clifford+T | 2.468e-2 | 4.616e-4 | 2,860 |

This is a variance-planning estimate, not a guarantee, and it ignores the
uncertainty in the estimated standard deviation. It supports an order-100
budget for the two noiseless cases. The unmodified noisy sampler would require
far more passes to match the deterministic error reliably, which quantifies the
need for noise-aware sampling rather than merely observing it qualitatively.

Raw measurements are in `optimal_pps_single_pass.csv` and
`optimal_pps_100_pass.csv`; the earlier 20-pass file is retained for exact
reproduction of the initial experiment.
