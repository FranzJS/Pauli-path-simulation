# Pauli-path simulation benchmark

The active benchmark now compares only:

- dense state-vector reference
- breadth-first Pauli propagation with a per-layer relative L1-mass cutoff of `1e-4`

Both retained circuit families use 21 qubits:

- random 1D Clifford+T brickwork: 12 layers, T density 0.70
- Trotterized nonintegrable Ising dynamics: 5 steps

At 21 qubits, the state vector contains `2^21` complex doubles and uses 32 MB, crossing the benchmark's 20 MB stopping threshold.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/pauli_benchmark 120
```

The argument is the L1-BFS timeout in seconds per circuit. Results are recorded in `results/benchmark_retained.csv`. Memory figures for BFS are representation-based frontier estimates; state-vector memory is the raw amplitude-array size.
