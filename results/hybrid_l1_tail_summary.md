# Hybrid L1-heavy-hitters + randomized-tail benchmark

## Algorithm

At each truncation event (every four reverse-propagated non-Clifford `RZ` gates), coefficients are sorted by absolute value. The smallest prefix whose total absolute mass is at most `cutoff * ||x||_1` is the tail; the remaining `K` coordinates are retained deterministically. The benchmark then retains exactly `ceil(alpha*K)` additional distinct tail coordinates for `alpha in {0.05, 0.10, 0.20}`.

Tail coordinates are chosen by fixed-size PPS sampling without replacement. Inclusion probabilities satisfy `pi_i = min(1, lambda*|x_i|)` and `sum_i pi_i = ceil(alpha*K)`; systematic PPS realizes the sample in linear time. The primary benchmark propagates selected coordinates with their original coefficients, matching the literal heuristic requested. An optional Horvitz–Thompson mode is implemented, but it is not used in the reported comparison because its reweighted tail can feed back into later heavy sets and destroy the intended resource profile.

Each randomized pass was executed in a fresh process. Memory is the repository's existing active-map estimate of 48 bytes per stored Pauli key, not operating-system RSS.

## Single-pass comparison

| Case | Method | Tail | Status | Runtime (s) | Runtime / BFS | Memory (MB) | Memory / BFS | Absolute error |
|---|---:|---:|---|---:|---:|---:|---:|---:|
| ct_n20_d12 | BFS L1 | — | ok | 6.419 | 1.00× | 37.77 | 1.00× | 7.917e-04 |
| ct_n20_d12 | Hybrid | 5% | ok | 11.535 | 1.80× | 49.28 | 1.30× | 4.054e-05 |
| ct_n20_d12 | Hybrid | 10% | ok | 11.951 | 1.86× | 50.28 | 1.33× | 2.442e-15 |
| ct_n20_d12 | Hybrid | 20% | ok | 15.933 | 2.48× | 54.65 | 1.45× | 2.442e-15 |
| ising_n20_s12 | BFS L1 | — | ok | 24.722 | 1.00× | 29.67 | 1.00× | 4.502e-04 |
| ising_n20_s12 | Hybrid | 5% | support_cap | 28.486 | 1.15× | 36.62 | 1.23× | — |
| ising_n20_s12 | Hybrid | 10% | support_cap | 27.972 | 1.13× | 36.62 | 1.23× | — |
| ising_n20_s12 | Hybrid | 20% | support_cap | 20.461 | 0.83× | 32.04 | 1.08× | — |
| ct_depol_n20_d12 | BFS L1 | — | ok | 0.067 | 1.00× | 0.87 | 1.00× | 4.616e-04 |
| ct_depol_n20_d12 | Hybrid | 5% | ok | 0.800 | 11.91× | 4.67 | 5.36× | 6.288e-05 |
| ct_depol_n20_d12 | Hybrid | 10% | ok | 1.639 | 24.40× | 9.59 | 11.02× | 3.141e-05 |
| ct_depol_n20_d12 | Hybrid | 20% | ok | 6.408 | 95.41× | 31.75 | 36.48× | 2.606e-05 |

The Ising hybrid rows are safety-stop diagnostics: 5% and 10% reached 800,001 active terms; 20% reached 700,001 active terms. No final estimate was produced, so no error is reported.

## Ten-pass cumulative means

| Case | Tail | Mean runtime/pass (s) | Total runtime (s) | Peak memory (MB) | 1-pass error | 10-pass mean error | Standard error at 10 passes |
|---|---:|---:|---:|---:|---:|---:|---:|
| ct_depol_n20_d12 | 5% | 0.695 | 6.953 | 5.03 | 6.288e-05 | 2.835e-05 | 7.527e-06 |
| ct_depol_n20_d12 | 10% | 1.444 | 14.438 | 9.73 | 3.141e-05 | 3.066e-05 | 8.415e-06 |
| ct_depol_n20_d12 | 20% | 5.593 | 55.933 | 32.43 | 2.606e-05 | 9.918e-06 | 4.423e-06 |
| ct_n20_d12 | 5% | 10.458 | 104.579 | 49.30 | 4.054e-05 | 5.332e-05 | 1.916e-05 |
| ct_n20_d12 | 10% | 10.767 | 107.672 | 50.28 | 2.442e-15 | 5.924e-06 | 2.418e-06 |
| ct_n20_d12 | 20% | 12.753 | 127.531 | 54.67 | 2.442e-15 | 2.054e-15 | 1.850e-17 |

No ten-pass Ising result is reported because the first pass did not complete within the support safety envelope.

## Main findings

1. **The single-pass runtimes are not generally "very similar".** On Clifford+T the hybrid costs about 1.8–2.5× BFS and 1.30–1.45× the estimated peak map memory. On the noisy case the BFS baseline is extremely small, while retained random-tail terms persist; the hybrid is about 12–95× slower and 5–36× larger.
2. **Accuracy improves dramatically on the Clifford+T instance.** A 5% tail reduces the selected single-pass error from `7.92e-4` to `4.05e-5`; 10% and 20% happen to recover the reference to numerical precision for the first seed. Across ten seeds, the 10% cumulative-mean error is `5.92e-6`, and 20% remains at numerical precision.
3. **The Ising circuit is the blocking case.** Even the 5% tail causes support growth beyond the chosen near-baseline safety envelope before the circuit finishes. This is consistent with repeated truncation/branching feedback: sampled tail coordinates can become part of later deterministic heavy sets.
4. **Averaging ten literal passes is not an unbiased Monte Carlo correction.** Selected tail coordinates keep their original coefficients, so averaging reduces random variation around a biased heuristic. The error is therefore not guaranteed to improve monotonically; for Clifford+T at 5%, the ten-pass error (`5.33e-5`) is slightly worse than the first-pass error (`4.05e-5`).
5. **The randomized tail can be useful, but it is circuit-dependent.** It is very effective on this Clifford+T observable, moderately useful on the noisy circuit, and operationally unstable on the Ising benchmark under this direct feedback design.

## Reproduction

Build the branch in Release mode, then use the one-pass executable from a fresh process for each seed. The committed raw CSV contains all 60 completed randomized passes and the three Ising safety-stop rows. The cumulative CSV records the running mean and error after every pass.
