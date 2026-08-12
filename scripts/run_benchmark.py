#!/usr/bin/env python3
"""Run each benchmark case in fresh processes and aggregate stable timings."""

from __future__ import annotations

import csv
import statistics
import subprocess
import sys
from pathlib import Path

RAW_COLUMNS = [
    "family",
    "case",
    "method",
    "truncation_strategy",
    "qubits",
    "l1_cutoff",
    "maximum_support",
    "minimum_magnitude",
    "estimate",
    "reference",
    "reference_method",
    "reference_id",
    "absolute_error",
    "runtime_s",
    "peak_support_terms",
    "peak_pre_truncation_terms",
    "peak_post_truncation_terms",
    "truncation_events",
    "peak_vector_capacity_terms",
    "peak_merge_entries",
    "peak_merge_capacity_slots",
]

OUTPUT_COLUMNS = [
    "family",
    "case",
    "method",
    "truncation_strategy",
    "qubits",
    "l1_cutoff",
    "maximum_support",
    "minimum_magnitude",
    "estimate",
    "reference",
    "reference_method",
    "reference_id",
    "absolute_error",
    "repetitions",
    "median_runtime_s",
    "min_runtime_s",
    "max_runtime_s",
    "peak_rss_mb",
    "peak_support_terms",
    "peak_pre_truncation_terms",
    "peak_post_truncation_terms",
    "truncation_events",
    "peak_vector_capacity_terms",
    "peak_merge_entries",
    "peak_merge_capacity_slots",
]

# Explicitly reproduce the three historical fixed benchmark configurations.
# The C++ executable is now general and accepts these parameters directly.
CASES = [
    (20, 12, "clifford_t", 0.00625),
    (20, 12, "ising", 0.00005),
    (20, 12, "clifford_t_depol", 0.034),
]


def run_once(
    executable: Path,
    configuration: tuple[int, int, str, float],
) -> dict[str, str]:
    qubits, layers, model, l1_cutoff = configuration
    completed = subprocess.run(
        [
            "/usr/bin/time",
            "-f",
            "__MAX_RSS_KB__=%M",
            str(executable),
            str(qubits),
            str(layers),
            model,
            "l1",
            str(l1_cutoff),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    values = next(csv.reader([completed.stdout.strip()]))
    if len(values) != len(RAW_COLUMNS):
        raise RuntimeError(f"unexpected benchmark output: {completed.stdout!r}")
    row = dict(zip(RAW_COLUMNS, values, strict=True))
    rss_line = next(
        line
        for line in completed.stderr.splitlines()
        if line.startswith("__MAX_RSS_KB__=")
    )
    row["peak_rss_mb"] = str(float(rss_line.split("=", 1)[1]) / 1024.0)
    return row


def main() -> int:
    build_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "build")
    repetitions = int(sys.argv[2] if len(sys.argv) > 2 else "5")
    output_path = Path(
        sys.argv[3] if len(sys.argv) > 3 else "results/benchmark_final.csv"
    )
    if repetitions <= 0:
        raise ValueError("repetitions must be positive")

    executable = build_dir / "pauli_benchmark"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []

    for case_index, configuration in enumerate(CASES):
        runs = [
            run_once(executable, configuration) for _ in range(repetitions)
        ]
        first = runs[0]
        runtimes = [float(run["runtime_s"]) for run in runs]
        rss_values = [float(run["peak_rss_mb"]) for run in runs]

        invariant_columns = [
            column for column in RAW_COLUMNS if column != "runtime_s"
        ]
        for run in runs[1:]:
            for column in invariant_columns:
                if run[column] != first[column]:
                    raise RuntimeError(
                        f"non-deterministic {column} in case {case_index}: "
                        f"{first[column]} != {run[column]}"
                    )

        rows.append(
            {
                **{column: first[column] for column in invariant_columns},
                "repetitions": repetitions,
                "median_runtime_s": f"{statistics.median(runtimes):.9f}",
                "min_runtime_s": f"{min(runtimes):.9f}",
                "max_runtime_s": f"{max(runtimes):.9f}",
                "peak_rss_mb": f"{max(rss_values):.3f}",
            }
        )

    with output_path.open("w", newline="") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
