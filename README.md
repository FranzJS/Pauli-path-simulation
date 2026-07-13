# Pauli-path simulation benchmarks

A modular C++20 benchmark comparing simulation methods on two entangling circuit families:

- random 1D Clifford+T brickwork: 14 qubits, 12 layers, T density 0.70
- Trotterized nonintegrable Ising dynamics: 16 qubits, 5 steps

Circuit definitions are separate from simulation methods.

## Methods

- dense state-vector reference
- exact sparse BFS propagation
- BFS with per-layer relative L1-mass truncation (discard at most 1e-4 of each layer's L1 coefficient mass)
- BFS with top-coefficient memory caps of 2, 4, and 8 MB
- exhaustive DFS path enumeration
- importance-sampled Monte Carlo paths

## Build and run

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/pauli_benchmark 25
```

Results are written to `results/benchmark_truncated.csv`. Memory figures are representation-based algorithmic estimates, not measured process RSS. Timed-out DFS partial sums are intentionally recorded as missing estimates.
