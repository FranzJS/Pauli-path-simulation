# Pauli-path simulation benchmarks

A modular C++20 benchmark comparing simulation methods on two entangling circuit families:

- random 1D Clifford+T brickwork: 14 qubits, 12 layers, T density 0.70
- Trotterized nonintegrable Ising dynamics: 16 qubits, 5 steps

Circuit definitions are separate from simulation methods in `src/v1_circuits.cpp` and `src/v1_methods.cpp`.

## Methods

- dense state-vector reference
- exact sparse Heisenberg Pauli propagation
- exhaustive depth-first Pauli-path enumeration
- importance-sampled Monte Carlo Pauli paths

## Build and run

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/pauli_benchmark 25
```

The argument is the time budget in seconds per method and case. The current benchmark records runtime, estimated algorithmic memory, and absolute error. Results are written to `results/benchmark_v1.csv`. DFS reaches the time cap on both pushed cases; exact sparse finishes with about 14 MB estimated frontier memory. The memory figures are representation-based estimates, not measured process RSS.
