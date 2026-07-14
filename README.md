# Pauli-path simulation benchmark

Stable benchmark v1 compares dense state-vector references and L1-truncated breadth-first Pauli propagation.

## Fixed cases

All cases use 20 qubits. L1 pruning occurs after every four non-Clifford `RZ` gates.

- Clifford+T: 12 brickwork layers, T density 0.70, seed 20260715, observable `X7 Z8 X9 X10`, cutoff 0.00625.
- Nonintegrable Ising: 12 Trotter steps, dt 0.12, J 1.0, hx 0.91, hz 0.37, observable `Z9 Z10`, cutoff 0.00005.
- Noisy Clifford+T: matched Clifford+T circuit with single-qubit depolarizing noise p=0.05 after each layer, observable `X7 Z8 X9 X10`, cutoff 0.034.

The cutoffs target absolute errors in the mid-10^-4 range.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/pauli_benchmark 120
```

Results are written to `results/benchmark_final.csv`. BFS memory is an algorithmic estimate based on 48 bytes per peak map term, not measured process RSS; transient allocator overhead may be higher.
