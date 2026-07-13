# Benchmark 0: scaling decomposition

This benchmark separates circuit definitions from simulation methods. The initial circuit family applies independent `H-Rz(theta)-H = Rx(theta)` rotations to `k` active qubits in an `n`-qubit register and measures `Z_0 ... Z_{k-1}` on `|0^n>`. The analytic reference is `cos(theta)^k`.

## Methods

- dense state-vector reference
- exact sparse Heisenberg Pauli propagation
- exhaustive depth-first path enumeration
- coefficient-magnitude importance-sampled Monte Carlo

Each method-case pair has a 25-second safety budget. The standard run also uses finite guards: 500,000 active terms for exact sparse propagation and one million Monte Carlo samples.

## Main findings

- With `k=12` fixed, increasing total width from 12 to 64 had negligible impact on the Pauli-path methods because the additional qubits were inactive.
- Exact sparse storage followed the designed `2^k` growth and reached the 500,000-term guard at `k=20`.
- DFS traversed the same `2^k` terminal paths but was faster on this collision-free family because it avoided hash-table and frontier overhead. This should not be generalized to circuits with strong path merging.
- Monte Carlo retained fixed memory, but at `theta=pi/4` its absolute error became comparable to the exponentially small signal. At `k=18`, one million samples produced error of order `10^-3` for a reference value of about `1.95e-3`.
- The angle sweep confirms that estimator variance, rather than register width alone, is the central Monte Carlo scaling hazard.

The canonical machine-readable result subset is `results/benchmark_summary.csv`. The TeX report contains the formal interpretation; the compiled PDF is provided with the benchmark deliverables.
