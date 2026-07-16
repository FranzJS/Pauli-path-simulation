# L1 heavy hitters plus randomized tail benchmark

## Experiment

This comparison uses the fixed 20-qubit benchmark from draft PR #2:

- `ct_n20_d12`, cutoff `0.00625`
- `ising_n20_s12`, cutoff `0.00005`
- `ct_depol_n20_d12`, cutoff `0.034`

The deterministic baseline is the repository's breadth-first Pauli propagation with merging and L1 truncation after every four non-Clifford `RZ` gates.

For a runtime-matched randomized pass, one deterministic calibration run records the number of retained heavy coordinates `K_t` at every truncation event. A randomized pass then:

1. keeps the largest `K_t` coefficients exactly;
2. samples exactly `ceil(alpha K_t)` distinct coordinates from the remaining tail, with fixed-size PPS probabilities proportional to coefficient magnitude;
3. propagates the selected tail coordinates with their original coefficients.

The tested tail ratios are `alpha = 0.05, 0.10, 0.20`. The method is deliberately **not** Horvitz-Thompson reweighted: reweighting caused sampled tail coefficients to become future heavy hitters and destroyed the requested runtime matching. Consequently, this is a biased randomized truncation heuristic. The multi-pass estimate is the arithmetic mean of independent passes.

Runtime for a randomized pass excludes the one-time deterministic calibration pass. For a one-off use, add the baseline runtime once. Memory is the repository's algorithmic estimate `48 bytes × peak pre-truncation support`, not process RSS.

## Single-pass comparison

| Case | Method | Runtime (s) | Memory (MB) | Absolute error | Peak pre-truncation terms |
|---|---:|---:|---:|---:|---:|
| ct_n20_d12 | BFS+L1 | 4.396 | 39.01 | 1.063e-03 | 852,188 |
| ct_n20_d12 | Heavy+tail 5% | 4.130 | 41.03 | 1.756e-03 | 896,228 |
| ct_n20_d12 | Heavy+tail 10% | 6.804 | 42.89 | 6.209e-04 | 936,912 |
| ct_n20_d12 | Heavy+tail 20% | 5.365 | 46.02 | 7.551e-04 | 1,005,398 |
| ising_n20_s12 | BFS+L1 | 13.623 | 29.67 | 4.421e-04 | 648,057 |
| ising_n20_s12 | Heavy+tail 5% | 19.379 | 31.36 | 3.861e-04 | 684,978 |
| ising_n20_s12 | Heavy+tail 10% | 19.779 | 32.90 | 3.071e-04 | 718,781 |
| ising_n20_s12 | Heavy+tail 20% | 20.234 | 35.97 | 2.669e-04 | 785,758 |
| ct_depol_n20_d12 | BFS+L1 | 0.062 | 0.87 | 4.616e-04 | 18,949 |
| ct_depol_n20_d12 | Heavy+tail 5% | 0.040 | 0.91 | 4.616e-04 | 19,918 |
| ct_depol_n20_d12 | Heavy+tail 10% | 0.041 | 0.96 | 1.102e-03 | 21,073 |
| ct_depol_n20_d12 | Heavy+tail 20% | 0.048 | 1.04 | 1.237e-03 | 22,822 |

## Progressive average after 10 passes

| Case | Tail | Mean runtime/pass (s) | Total runtime (s) | Max memory (MB) | Error of 10-pass mean | Empirical SD | Empirical SE |
|---|---:|---:|---:|---:|---:|---:|---:|
| ct_n20_d12 | 5% | 4.574 | 45.740 | 41.40 | 5.450e-04 | 8.130e-04 | 2.571e-04 |
| ct_n20_d12 | 10% | 5.036 | 50.365 | 42.92 | 1.967e-04 | 5.158e-04 | 1.631e-04 |
| ct_n20_d12 | 20% | 6.063 | 60.633 | 46.02 | 3.357e-04 | 3.365e-04 | 1.064e-04 |
| ising_n20_s12 | 5% | 17.999 | 179.993 | 31.39 | 4.002e-04 | 1.595e-05 | 5.042e-06 |
| ising_n20_s12 | 10% | 18.921 | 189.211 | 32.98 | 3.309e-04 | 2.420e-05 | 7.654e-06 |
| ising_n20_s12 | 20% | 19.741 | 197.410 | 35.99 | 2.830e-04 | 1.729e-05 | 5.467e-06 |
| ct_depol_n20_d12 | 5% | 0.037 | 0.367 | 0.91 | 4.410e-04 | 6.334e-04 | 2.003e-04 |
| ct_depol_n20_d12 | 10% | 0.040 | 0.396 | 0.97 | 7.381e-04 | 6.802e-04 | 2.151e-04 |
| ct_depol_n20_d12 | 20% | 0.043 | 0.432 | 1.05 | 2.048e-06 | 7.886e-04 | 2.494e-04 |

## Main findings

- **Clifford+T:** deterministic BFS error was `1.063e-03`. The 10% tail gave the best 10-pass result, `1.967e-04`, at `5.04` seconds per pass and about 10% more peak frontier memory. The prefix error was non-monotone and reached `1.781e-05` after 5 passes before rising again.
- **Ising:** all three randomized single passes improved on BFS. The 20% tail was best, with single-pass error `2.669e-04` and 10-pass error `2.830e-04` versus BFS `4.421e-04`. Averaging did not materially reduce the error, indicating a stable truncation bias.
- **Noisy Clifford+T:** results were seed-sensitive. The 20% tail's 10-pass mean happened to land within `2.048e-06` of the reference, but its empirical per-pass SD was `7.886e-04`. This near-cancellation should not be treated as robust evidence from only 10 passes.
- The calibrated support cap worked: peak frontier memory rose roughly with `1 + alpha`, while pass runtime stayed in the same single-core regime. Ising incurred the largest overhead, about 1.4–1.5× BFS.

## Interpretation

The method is promising as a low-complexity accuracy enhancement, especially for the Clifford+T and Ising cases. It does not behave like an unbiased Monte Carlo estimator: error trajectories can plateau or move away from the reference as more passes are averaged. A larger seed study is needed before selecting a universal tail ratio. On these three cases, 10% is the strongest default for Clifford+T, while 20% performs best for Ising; the noisy result is too variable for a firm choice.

The repository's L1 truncation comparator only orders by coefficient magnitude. Equal-magnitude boundary terms therefore depend on hash-map iteration order, as already noted in the repository audit. The baseline and randomized runs here were produced by the same executable, compiler, and machine, so their comparison is internally consistent even though the absolute baseline values can differ slightly from previously committed CSVs.

## Environment

- Linux x86_64
- GCC 14.2.0
- Release-style compile: `-O3 -march=native`
- Single-threaded runs
