# Pauli-path simulation benchmarks

A modular C++20 benchmark for comparing circuit simulation methods. Circuit families are defined independently in `src/circuits.cpp`; simulation methods are implemented in `src/methods.cpp`.

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

The argument is the time budget in seconds **per method and case**. Results are written to `results/benchmark.csv`. The formal report source is in `report/benchmark_report.tex`; `report/summary.md` provides a browsable repository summary. The compiled PDF can be regenerated with `pdflatex`.
