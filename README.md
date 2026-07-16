# Pauli-path simulation benchmark

This repository contains one reproducible benchmark for **breadth-first Heisenberg Pauli propagation with exact path merging and periodic L1-mass truncation**.

The implementation deliberately has one canonical Pauli propagation kernel. Clifford gates and diagonal noise are applied in place to a unique frontier vector. Only branching `RZ` gates allocate a merge table for exact collision handling. Trigonometric factors are precomputed once per gate, and the merge table is a reusable flat open-addressing structure specialized for Pauli keys.

## Fixed benchmark

All cases use 20 qubits. Truncation occurs after every four processed `RZ` gates. At each truncation, the smallest coefficients are removed while their cumulative absolute mass does not exceed `cutoff * current_L1_mass`.

| Case | Parameters | Observable | L1 cutoff | Reference |
|---|---|---|---:|---:|
| Clifford+T | 12 brickwork layers, T density 0.70, seed 20260715 | `X7 Z8 X9 X10` | 0.00625 | 0.329520762689497 |
| Nonintegrable Ising | 12 Trotter steps, dt 0.12, J 1.0, hx 0.91, hz 0.37 | `Z9 Z10` | 0.00005 | 0.714108793456576 |
| Noisy Clifford+T | matched circuit, single-qubit depolarizing p=0.05 after every layer | `X7 Z8 X9 X10` | 0.034 | 0.0663217307060528 |

The noiseless references are dense statevector results. The noisy reference is a converged Pauli calculation.

## Build and test

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

For a machine-specific benchmark build:

```bash
cmake -S . -B build-native \
  -DCMAKE_BUILD_TYPE=Release \
  -DPAULI_NATIVE_OPT=ON
cmake --build build-native --parallel
```

## Reproduce the benchmark

Each case is run in a separate process so that peak RSS is not contaminated by earlier cases.

```bash
./scripts/run_benchmark.py build-native 5 results/benchmark_final.csv
```

The committed result file records the median and minimum of five repetitions, peak process RSS, and logical frontier diagnostics. Runtime and RSS are machine-dependent; estimates and errors should reproduce up to floating-point and tied-truncation ordering effects.

## Committed benchmark result

The committed CSV was produced on a 56-core x86_64 virtual machine using GCC 14.2, CMake 3.31.6, `Release`, and `PAULI_NATIVE_OPT=ON`. Each timing is the median of five fresh-process runs.

| Case | Median runtime | Peak RSS | Absolute error |
|---|---:|---:|---:|
| Clifford+T | 0.861 s | 71.9 MB | 6.91e-4 |
| Nonintegrable Ising | 3.292 s | 69.3 MB | 4.50e-4 |
| Noisy Clifford+T | 0.0136 s | 9.59 MB | 4.62e-4 |

The canonical kernel is substantially faster than the previous node-based implementation while preserving the benchmark definition. The exact retained support can differ slightly from older runs because truncation ties are now broken deterministically by the Pauli key.
