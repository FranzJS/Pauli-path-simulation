# Resampling BFS implementation audit

## Main finding

The current implementation uses the same integer `K` for three different roles:

1. the number of particle copies drawn or rounded;
2. the trigger `current.size() > K` for compression;
3. the nominal retained-memory estimate `48*K`.

These quantities are not equal after repeated particle copies of the same Pauli are merged. `rebuild` stores one map entry per nonzero count, not one entry per particle. Therefore the previous `nominal_mb` values are particle budgets, not retained-map memory budgets, and the implementation does not literally retain `K` distinct Pauli coefficients.

The calibrated CSV did use observed pre-resampling peak support, so its peak-memory comparisons are not entirely invalid. However, the stochastic maps after resampling were often much smaller than the deterministic BFS maps after L1 truncation. The earlier interpretation that the algorithms tracked approximately the same number of heavy coefficients was incorrect.

## Correctness checks

- The H, S, CNOT, RZ and depolarizing transformations agree with the deterministic BFS implementation.
- Ordinary multinomial, residual multinomial and pairwise dependent rounding preserve the required conditional means.
- A 200,000-trial synthetic test of the dependent-rounding primitive reproduced all target marginals within Monte Carlo error and preserved the requested total population exactly.
- Full-circuit sample means were consistent with the references within their empirical standard errors.
- Repeated copies of one Pauli are collapsed into one coefficient before further propagation; they are not propagated separately inside a replica.

No algebraic bias bug was found.

## Representative diagnostics

On one audited seed, residual dependent rounding gave approximately:

- Clifford+T, nominal 16 MB particles: 32.75 MB pre-resampling peak, 13.50 MB maximum post-resampling map, four stochastic compression events. At the last event, 46,389 of 715,381 coordinates had expected count at least one; these coordinates carried about 29% of the L1 mass.
- Ising, nominal 10 MB particles: 33.81 MB pre-resampling peak, 9.70 MB maximum post-resampling map, 21 stochastic compression events. At a late event, 9,857 of 507,708 coordinates had expected count at least one and carried about 25% of the L1 mass.
- Noisy Clifford+T, nominal 2.8 MB particles: 7.28 MB pre-resampling peak, 2.56 MB maximum post-resampling map, six stochastic compression events.

For comparison, an instrumented L1-BFS run retained maximum post-truncation maps of roughly 34.6 MB for Clifford+T, 24.4 MB for Ising and 0.68 MB for the noisy case. Exact values can vary slightly at tied cutoff boundaries because the hash-map iteration order affects which equal-magnitude coefficient is removed last.

## Avoidable repeated work

- `advance` reserves twice the current support for every gate, although only RZ gates can branch. This can be reduced for H, S, CNOT and DEPOL.
- Ordinary and residual multinomial currently make one `discrete_distribution` draw per particle. An alias-table or grouped count sampler can reduce sampling overhead. This is sampling overhead, not duplicate Pauli propagation.
- Every replica recomputes the identical deterministic prefix before the first stochastic compression. It can be computed once and distributed or forked with copy-on-write. This reduces aggregate core-work but normally not one-wave wall time.
- After the first stochastic event, sharing overlapping heavy-key propagation across replicas would require communication and cross-replica merging, sacrificing the embarrassingly parallel design.

## Required redesign

A memory-faithful implementation must separate:

- `B`: maximum number of distinct stored Pauli keys / actual support-memory cap;
- `P`: particle resolution used to control coefficient quantization variance.

Compression should be triggered by support relative to `B`, not by `P`. A resampler must then produce at most `B` distinct keys. The most direct exact-cap construction is fixed-size PPS / dependent inclusion rounding with marginal inclusion probabilities and Horvitz-Thompson coefficients. The present particle-count residual rounding does not by itself enforce a distinct-key cap when `P != B`.

The previous 100-core projections should therefore be treated as results for the `B=P` prototype, not as a definitive test of the intended memory-capped heavy-coefficient algorithm.
