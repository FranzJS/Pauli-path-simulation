# Pauli-path simulation benchmark

The active benchmark compares only:

- dense state-vector reference
- breadth-first Pauli propagation with a per-layer relative L1-mass cutoff

Both circuit families use 21 qubits:

- random 1D Clifford+T brickwork: 12 layers, T density 0.70
- Trotterized nonintegrable Ising dynamics: 12 steps

The L1 cutoff was doubled from `1e-4` until the full Clifford+T run completed within 120 seconds. The first successful tested value was `0.0128`; all tested values through `0.0064` exceeded the time window. The same global cutoff is used for both circuits.

At 21 qubits, the state vector contains `2^21` complex doubles and uses 32 MB.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/pauli_benchmark 120
```

Results are recorded in `results/benchmark_retained.csv`. BFS memory is a retained-frontier estimate; state-vector memory is the raw amplitude-array size.
