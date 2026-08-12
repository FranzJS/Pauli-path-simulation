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

## Reference registry

Exact results and their identifying metadata live in
[`references/reference_registry.json`](references/reference_registry.json), not
in the C++ circuit source. A record is selected only when all configuration
fields match: model, qubits, layers/Trotter steps, circuit-generation version,
seed, physical parameters, and both observable masks. Thus a different size or
seed cannot accidentally inherit an unrelated reference.

Validate the registry after adding or editing a record:

```bash
python3 scripts/validate_reference_registry.py
```

The validator checks schema version 2, required fields and parameter ranges,
model-specific fields, Pauli masks, finite results, unique IDs, and duplicate or
conflicting configurations. CTest runs this validation automatically when a
Python interpreter is available.

By default, C++ programs use the registry in this source tree. To run against a
different registry without rebuilding, set its path for that process:

```bash
PAULI_REFERENCE_REGISTRY=/path/to/references.json \
  ./build/pauli_benchmark 20 12 clifford_t l1 0.00625
```

The public lookup interface is declared in
[`reference_registry.hpp`](include/pauli_bench/reference_registry.hpp). A valid
registry with no exact match produces a `nan` value/error, `unavailable`
method, and empty reference ID; an unreadable, malformed, or ambiguous registry
stops the run with an error.

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
| `clifford_t_identity` | optional shallow W, then Clifford+T U and its exact inverse | same noiseless Clifford+T parameters |
| `ising_identity` | optional shallow W, then Trotterized Ising U and its exact inverse | same Ising parameters |

`LAYERS` means brickwork layers for the Clifford models and Trotter steps for
Ising. For an `_identity` model it is the forward depth; the constructor appends
the same number of inverse logical layers, giving `2 * LAYERS` total.
`L1_CUTOFF` is the maximum fraction of the current frontier L1 mass that may be
removed at each truncation event; it must be between 0 and 1. It is a mass-error
budget, not a direct maximum number of retained terms.

With `support`, every term satisfying `abs(coefficient) < MIN_MAGNITUDE` is
removed first. If more than `MAX_SUPPORT` terms remain, only the
`MAX_SUPPORT` largest-magnitude terms are retained. Equality with the threshold
is retained unless the support cap removes it.

The Clifford models support 4--64 qubits and Ising supports 2--64. The
observable is centered as the width changes. Only the original 20-qubit,
12-layer/step standard circuits have stored registry values; other standard
configurations print `nan,unavailable,,nan` for reference, reference method,
empty reference ID, and absolute error. Prefixed identity circuits instead use
the exact lightcone-statevector expectation of the shallow `W` circuit. CSV
output includes the reference ID so generated results retain the provenance
link.

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
  --passes 100 --batches 5 \
  --output results/support_pps_convergence_5x100 --workers 8

./scripts/run_pps_convergence.py build-native l1 \
  0.00625,0.00005,0.034 \
  --passes 100 --output results/l1_pps_convergence_100 --workers 8
```

Each worker launches a separate single-threaded C++ process. Results are reduced
in batch/pass-index order, so a fixed configuration and seed sequence produce
the same estimates and running statistics regardless of worker count; timing
and peak-RSS fields remain machine- and load-dependent. `total_runtime_s` and
`cumulative_pps_runtime_s` are sums of per-pass compute times, not elapsed wall
time. Choose a smaller worker count when simultaneous frontier memory would
exceed available RAM or when memory-bandwidth contention stops improving wall
time.

With `--batches B`, the convergence SVG draws each batch's running-mean
absolute-error curve transparently and the pointwise arithmetic average of
those error curves in full color. It does not draw a standard-error curve; the
deterministic BFS error remains a dashed horizontal line. The former
`run_support_pps_convergence.py` path remains as a compatibility wrapper.
The pass-axis tick spacing is selected automatically from 1, 2, or 5 times a
power of ten, keeping approximately six labels readable even for runs with
hundreds or thousands of passes.

## Identity-circuit convergence

[`run_identity_pps_convergence.py`](scripts/run_identity_pps_convergence.py)
runs the same convergence experiment at a chosen width for the two noiseless
families. For `n` qubits it constructs a shallow circuit `W`, then `n` ordinary
`U` layers or Trotter steps, followed by the exact gatewise inverse of `U`.
Thus the full circuit has depth `W_DEPTH + 2n` and acts exactly like `W`.
`W_DEPTH` defaults to 5.

The exact reference is
`<0|W^dagger O W|0>`. On first use of a new width, W depth, or Clifford circuit
seed, `pauli_identity_reference` evaluates this on a reduced statevector
containing only the finite-depth backward lightcone. The script atomically adds
the result and its complete configuration to the JSON reference registry.
Later runs reuse that record.

Support-derived schedules:

```bash
python3 scripts/run_identity_pps_convergence.py \
  build 12 support 20000 1e-8 \
  --w-depth 5 --circuit-seed 20260715 \
  --passes 100 --batches 5 --workers 8 \
  --output results/identity_support_n12_5x100
```

L1-derived schedules, with separate Clifford+T and Ising cutoffs:

```bash
python3 scripts/run_identity_pps_convergence.py \
  build 12 l1 0.00625,0.00005 \
  --w-depth 5 --circuit-seed 20260715 \
  --passes 100 --batches 5 --workers 8 \
  --output results/identity_l1_n12_5x100
```

Both `--w-depth` and `--circuit-seed` may be omitted to use their defaults.
For Clifford+T, the argument is a master seed: it selects `U` directly and a
stable derived seed selects an independent `W` stream. This preserves the same
`U` realization as the corresponding unprefixed run without duplicating its
first layers in `W`. The Ising circuit is deterministic. `--registry PATH`
selects an alternative schema-version-2 registry and is useful when generated
references should not modify the repository registry.

The output prefix produces a CSV and a two-panel SVG. As in the general
convergence runner, transparent curves are individual batches, the solid curve
is their pointwise average running-mean absolute error, and the dashed line is
the deterministic BFS error.

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
