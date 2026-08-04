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

## Run a configurable simulation

The deterministic executable accepts the number of qubits, number of
layers/Trotter steps, model, and one of three K-determination strategies:

```bash
./build/pauli_benchmark QUBITS LAYERS MODEL l1 L1_CUTOFF
./build/pauli_benchmark QUBITS LAYERS MODEL \
  support MAX_SUPPORT MIN_MAGNITUDE
./build/pauli_benchmark QUBITS LAYERS MODEL schedule K0:K1:K2:...
```

For example, to run increasing Clifford+T sizes:

```bash
./build/pauli_benchmark 8  6 clifford_t l1 0.00625
./build/pauli_benchmark 12 6 clifford_t support 10000 1e-10
./build/pauli_benchmark 16 6 clifford_t support 100000 1e-12
```

Supported models are:

| Model | Meaning | Fixed model-specific parameters |
|---|---|---|
| `clifford_t` | noiseless brickwork | T density 0.70, seed 20260715 |
| `clifford_t_depol` | noisy brickwork | same circuit, depolarizing p=0.05 per qubit after each layer |
| `ising` | nonintegrable Ising | dt 0.12, J 1.0, hx 0.91, hz 0.37 |

`LAYERS` means brickwork layers for the Clifford models and Trotter steps for
Ising. `L1_CUTOFF` is the maximum fraction of the current frontier L1 mass that
may be removed at each truncation event; it must be between 0 and 1. It is a
mass-error budget, not a direct maximum number of retained terms.

With `support`, every term satisfying `abs(coefficient) < MIN_MAGNITUDE` is
removed first. If more than `MAX_SUPPORT` terms remain, only the
`MAX_SUPPORT` largest-magnitude terms are retained. Equality with the threshold
is retained unless the support cap removes it.

The Clifford models support 4--64 qubits and Ising supports 2--64. The
observable is centered as the width changes. Only the original 20-qubit,
12-layer/step circuits have stored reference values; other configurations print
`nan,unavailable,nan` for reference, reference method, and absolute error.

The configurable PPS interface is:

```bash
./build/pauli_magnitude_pps_benchmark \
  QUBITS LAYERS MODEL l1 L1_CUTOFF bfs

./build/pauli_magnitude_pps_benchmark \
  QUBITS LAYERS MODEL support MAX_SUPPORT MIN_MAGNITUDE schedule

./build/pauli_magnitude_pps_benchmark \
  QUBITS LAYERS MODEL schedule K0:K1:K2:... pps_ht SEED

./build/pauli_magnitude_pps_benchmark \
  QUBITS LAYERS MODEL support MAX_SUPPORT MIN_MAGNITUDE pps_ht SEED
```

K determination and K-term selection are separate. Deterministic propagation
supports `l1`, `support`, and `schedule`, then retains the K largest terms. PPS
supports `support` and `schedule`, then applies whole-frontier PPS/HT sampling.
`l1` is deliberately rejected for PPS. With `support`, PPS computes
`K = min(MAX_SUPPORT, count(abs(coefficient) >= MIN_MAGNITUDE))` from its own
current frontier at each event.

For one randomized pass, first obtain the deterministic support schedule:

```bash
schedule=$(./build/pauli_magnitude_pps_benchmark \
  12 6 ising support 10000 1e-10 schedule)

./build/pauli_magnitude_pps_benchmark \
  12 6 ising schedule "$schedule" pps_ht 12345

# Or let this PPS pass choose K from its own frontier at every event:
./build/pauli_magnitude_pps_benchmark \
  12 6 ising support 10000 1e-10 pps_ht 12345
```

## Reproduce the benchmark

Each case is run in a separate process so that peak RSS is not contaminated by earlier cases.

```bash
./scripts/run_benchmark.py build-native 5 results/benchmark_final.csv
```

The Python runner still reproduces the three fixed cases in the table above; it
now passes their qubits, layers, model, and cutoff explicitly to the configurable
C++ executable.

The committed result file records the median and minimum of five repetitions, peak process RSS, and logical frontier diagnostics. Runtime and RSS are machine-dependent; estimates and errors should reproduce up to floating-point and tied-truncation ordering effects.

## Committed benchmark result

The committed CSV was produced on a 56-core x86_64 virtual machine using GCC 14.2, CMake 3.31.6, `Release`, and `PAULI_NATIVE_OPT=ON`. Each timing is the median of five fresh-process runs.

| Case | Median runtime | Peak RSS | Absolute error |
|---|---:|---:|---:|
| Clifford+T | 0.861 s | 71.9 MB | 6.91e-4 |
| Nonintegrable Ising | 3.292 s | 69.3 MB | 4.50e-4 |
| Noisy Clifford+T | 0.0136 s | 9.59 MB | 4.62e-4 |

The canonical kernel is substantially faster than the previous node-based implementation while preserving the benchmark definition. The exact retained support can differ slightly from older runs because truncation ties are now broken deterministically by the Pauli key.

## Magnitude PPS with Horvitz--Thompson correction

The second simulator first runs deterministic L1 truncation once to record the
retained support budget `K` at every truncation event. An independent randomized
pass then samples exactly `K` coordinates from the entire frontier using

```text
pi_i = min(1, |c_i| / tau),  with  sum_i pi_i = K.
```

Coordinates with `pi_i = 1` are deterministic heavy hitters. The remaining
coordinates are drawn by fixed-size systematic probability-proportional-to-size
sampling and selected coefficients receive the Horvitz--Thompson factor
`1 / pi_i`. Thus every truncation, and therefore every complete pass, is
unbiased while matching the deterministic BFS support schedule.

Reproduce the single-pass and progressive comparisons with:

```bash
./scripts/run_optimal_pps_benchmark.py \
  build-native 5 \
  results/optimal_pps_single_pass.csv \
  results/optimal_pps_100_pass.csv \
  100 8
```

The final argument is the maximum number of concurrent PPS passes. If omitted,
it defaults to the number of logical CPUs reported by the operating system.
Timing repetitions used for the single-pass comparison remain sequential so
that their measurements are not distorted by competing workers. The general
convergence runner accepts either a support-derived or L1-derived deterministic
schedule and can execute PPS passes concurrently:

```bash
./scripts/run_pps_convergence.py build-native support 100000 1e-20 \
  --passes 100 --output results/support_pps_convergence_100 --workers 8

./scripts/run_pps_convergence.py build-native l1 \
  0.00625,0.00005,0.034 \
  --passes 100 --output results/l1_pps_convergence_100 --workers 8
```

Each worker launches a separate single-threaded C++ process. Results are reduced
in pass-index order, so a fixed configuration and seed sequence produce the
same estimates and running statistics regardless of worker count; timing and
peak-RSS fields remain machine- and load-dependent. `total_runtime_s` and
`cumulative_pps_runtime_s` are sums of per-pass compute times, not elapsed wall
time. Choose a smaller worker count when simultaneous frontier memory would
exceed available RAM or when memory-bandwidth contention stops improving wall
time.

The convergence SVG plots the PPS empirical standard error
`stdev(estimates) / sqrt(R)` on the same log scale as the running-mean absolute
error and BFS error. To avoid emphasizing unstable few-sample estimates, the SE
curve starts at `max(2, ceil(sqrt(total passes)))` (for example, pass 10 in a
100-pass run).
The pass-axis tick spacing is selected automatically from 1, 2, or 5 times a
power of ten, keeping approximately six labels readable even for runs with
hundreds or thousands of passes.

The committed outputs include the
[single-pass comparison](results/optimal_pps_single_pass.csv),
[20-pass trajectories](results/optimal_pps_20_pass.csv),
[100-pass trajectories](results/optimal_pps_100_pass.csv), and a concise
[benchmark report](results/optimal_pps_report.md). Randomized passes are
independent and can be executed in parallel; their estimates are averaged only
after all passes complete.

## Small Clifford+T PPS counterexample

The following 10-qubit example shows that the *realized* error of a finite PPS
mean need not decrease monotonically with the support budget, even when the
empirical standard error does:

```bash
cmake --build build --target \
  pauli_frontier_spectrum pauli_magnitude_pps_benchmark
python3 scripts/run_small_clifford_pps_example.py build --workers 8
```

With the default fixed circuit and 80 fixed PPS seeds, increasing `K` from
10,000 to 30,000 changes the final running-mean error from about 0.083 to 0.179,
while the empirical standard error falls from about 0.98 to 0.26. This is a
finite-sample reversal, not evidence that the larger-`K` estimator has larger
expected RMSE. The accompanying SVG plots both trajectories and the exact
pre-truncation Clifford+T coefficient staircase; vertical lines show where the
two support budgets cut the ranked spectrum.
