#!/usr/bin/env python3
"""Benchmark whole-frontier PPS variants with HT correction."""

from __future__ import annotations

import csv
import math
import os
import statistics
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

RAW_COLUMNS = [
    "family", "case", "method", "truncation_strategy", "tail_ratio",
    "seed", "qubits", "l1_cutoff", "maximum_support",
    "minimum_magnitude", "estimate", "reference", "reference_method",
    "absolute_error", "runtime_s", "peak_support_terms",
    "peak_pre_truncation_terms", "peak_post_truncation_terms",
    "truncation_events", "peak_vector_capacity_terms",
    "peak_merge_entries", "peak_merge_capacity_slots", "max_heavy_terms",
    "max_sampled_tail_terms", "total_sampled_tail_terms",
    "max_importance_multiplier", "max_post_truncation_abs_coefficient",
]

SINGLE_PASS_COLUMNS = [
    *[column for column in RAW_COLUMNS if column != "runtime_s"],
    "timing_repetitions", "median_runtime_s", "min_runtime_s",
    "max_runtime_s", "peak_rss_mb", "runtime_vs_bfs", "rss_vs_bfs",
]

PROGRESSIVE_COLUMNS = [
    "family", "case", "tail_ratio", "passes", "seed", "pass_estimate",
    "pass_absolute_error", "mean_estimate", "mean_absolute_error",
    "empirical_sd", "empirical_stderr", "pass_runtime_s",
    "mean_runtime_s", "total_runtime_s", "max_peak_rss_mb",
    "pass_max_importance_multiplier", "max_importance_multiplier",
    "pass_max_post_truncation_abs_coefficient",
    "max_post_truncation_abs_coefficient", "reference",
]

# Explicitly reproduce the three historical fixed benchmark configurations.
CASES = [
    (20, 12, "clifford_t", 0.00625),
    (20, 12, "ising", 0.00005),
    (20, 12, "clifford_t_depol", 0.034),
]


def run_once(
    executable: Path,
    configuration: tuple[int, int, str, float],
    method: str,
    seed: int = 0,
    retained_schedule: tuple[int, ...] = (),
) -> dict[str, str]:
    qubits, layers, model, l1_cutoff = configuration
    if method == "pps_ht":
        if not retained_schedule:
            raise ValueError("PPS runs require a retained-support schedule")
        command = [
            str(executable), str(qubits), str(layers), model,
            "schedule",
            ":".join(str(value) for value in retained_schedule),
            method,
            str(seed),
        ]
    else:
        command = [
            str(executable), str(qubits), str(layers), model,
            "l1", str(l1_cutoff), method,
        ]
    wrapper = Path(__file__).with_name("measure_process.py")
    completed = subprocess.run(
        [sys.executable, str(wrapper), *command],
        check=True,
        capture_output=True,
        text=True,
    )
    values = next(csv.reader([completed.stdout.strip()]))
    if len(values) != len(RAW_COLUMNS):
        raise RuntimeError(f"unexpected benchmark output: {completed.stdout!r}")
    row = dict(zip(RAW_COLUMNS, values, strict=True))
    rss_line = next(
        line for line in completed.stderr.splitlines()
        if line.startswith("__MAX_RSS_KB__=")
    )
    row["peak_rss_mb"] = str(float(rss_line.split("=", 1)[1]) / 1024.0)
    return row


def repeated_timing(
    executable: Path,
    configuration: tuple[int, int, str, float],
    method: str,
    repetitions: int,
    seed: int = 0,
    retained_schedule: tuple[int, ...] = (),
) -> tuple[dict[str, str], list[float], list[float]]:
    runs = [
        run_once(executable, configuration, method, seed, retained_schedule)
        for _ in range(repetitions)
    ]
    first = runs[0]
    invariant_columns = [column for column in RAW_COLUMNS if column != "runtime_s"]
    for run in runs[1:]:
        for column in invariant_columns:
            if run[column] != first[column]:
                raise RuntimeError(
                    f"non-deterministic {column}: {first[column]} != {run[column]}"
                )
    return (
        first,
        [float(run["runtime_s"]) for run in runs],
        [float(run["peak_rss_mb"]) for run in runs],
    )


def seed_for(case_index: int, pass_index: int) -> int:
    return 202707160300 + 1009 * case_index + pass_index


def default_workers() -> int:
    return os.cpu_count() or 1


def calibrate_schedule(
    executable: Path,
    configuration: tuple[int, int, str, float],
) -> tuple[int, ...]:
    qubits, layers, model, l1_cutoff = configuration
    completed = subprocess.run(
        [
            str(executable), str(qubits), str(layers), model,
            "l1", str(l1_cutoff), "schedule",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return tuple(int(value) for value in completed.stdout.strip().split(":"))


def single_pass_row(
    run: dict[str, str],
    runtimes: list[float],
    rss_values: list[float],
    repetitions: int,
    baseline_runtime: float,
    baseline_rss: float,
) -> dict[str, object]:
    median_runtime = statistics.median(runtimes)
    peak_rss = max(rss_values)
    return {
        **{
            column: run[column]
            for column in RAW_COLUMNS
            if column != "runtime_s"
        },
        "timing_repetitions": repetitions,
        "median_runtime_s": f"{median_runtime:.9f}",
        "min_runtime_s": f"{min(runtimes):.9f}",
        "max_runtime_s": f"{max(runtimes):.9f}",
        "peak_rss_mb": f"{peak_rss:.3f}",
        "runtime_vs_bfs": f"{median_runtime / baseline_runtime:.6f}",
        "rss_vs_bfs": f"{peak_rss / baseline_rss:.6f}",
    }


def main() -> int:
    build_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "build")
    repetitions = int(sys.argv[2] if len(sys.argv) > 2 else "5")
    single_path = Path(
        sys.argv[3]
        if len(sys.argv) > 3
        else "results/optimal_pps_single_pass.csv"
    )
    progressive_path = Path(
        sys.argv[4]
        if len(sys.argv) > 4
        else "results/optimal_pps_20_pass.csv"
    )
    pass_spec = sys.argv[5] if len(sys.argv) > 5 else "20"
    workers = int(sys.argv[6] if len(sys.argv) > 6 else default_workers())
    if "," in pass_spec:
        pass_counts = tuple(int(value) for value in pass_spec.split(","))
        if len(pass_counts) != 3:
            raise ValueError("per-case pass counts must contain three values")
    else:
        pass_counts = (int(pass_spec),) * 3
    method = "pps_ht"
    if repetitions <= 0:
        raise ValueError("timing repetitions must be positive")
    if any(pass_count <= 0 for pass_count in pass_counts):
        raise ValueError("pass counts must be positive")
    if workers <= 0:
        raise ValueError("workers must be positive")

    executable = build_dir / "pauli_magnitude_pps_benchmark"
    single_rows: list[dict[str, object]] = []
    progressive_rows: list[dict[str, object]] = []

    for case_index, configuration in enumerate(CASES):
        max_passes = pass_counts[case_index]
        retained_schedule = calibrate_schedule(executable, configuration)
        baseline, baseline_times, baseline_rss_values = repeated_timing(
            executable,
            configuration,
            "bfs",
            repetitions,
        )
        baseline_median = statistics.median(baseline_times)
        baseline_peak_rss = max(baseline_rss_values)
        single_rows.append(
            single_pass_row(
                baseline,
                baseline_times,
                baseline_rss_values,
                repetitions,
                baseline_median,
                baseline_peak_rss,
            )
        )

        first_seed = seed_for(case_index, 1)
        first, first_times, first_rss_values = repeated_timing(
            executable,
            configuration,
            method,
            repetitions,
            seed=first_seed,
            retained_schedule=retained_schedule,
        )
        first_median = statistics.median(first_times)
        first_peak_rss = max(first_rss_values)
        single_rows.append(
            single_pass_row(
                first,
                first_times,
                first_rss_values,
                repetitions,
                baseline_median,
                baseline_peak_rss,
            )
        )
        print(
            f"case {case_index}: BFS {baseline_median:.3f}s/"
            f"{baseline_peak_rss:.1f} MB; PPS {first_median:.3f}s/"
            f"{first_peak_rss:.1f} MB",
            flush=True,
        )

        estimates: list[float] = []
        runtimes: list[float] = []
        peak_rss_values: list[float] = []
        max_multiplier = 1.0
        max_abs_coefficient = 0.0
        reference = float(first["reference"])

        def run_pass(pass_index: int) -> dict[str, str]:
            return run_once(
                executable,
                configuration,
                method,
                seed=seed_for(case_index, pass_index),
                retained_schedule=retained_schedule,
            )

        pass_indices = range(2, max_passes + 1)
        worker_count = min(workers, max(1, max_passes - 1))
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            remaining_runs = iter(executor.map(run_pass, pass_indices))
            ordered_runs = [first, *remaining_runs]

        for pass_index, run in enumerate(ordered_runs, start=1):
            if pass_index == 1:
                runtime = first_median
                peak_rss = first_peak_rss
            else:
                runtime = float(run["runtime_s"])
                peak_rss = float(run["peak_rss_mb"])

            estimate = float(run["estimate"])
            if not math.isfinite(estimate):
                raise RuntimeError(
                    f"non-finite estimate for case {case_index}, "
                    f"pass {pass_index}"
                )
            estimates.append(estimate)
            runtimes.append(runtime)
            peak_rss_values.append(peak_rss)
            pass_multiplier = float(run["max_importance_multiplier"])
            pass_max_coefficient = float(
                run["max_post_truncation_abs_coefficient"]
            )
            max_multiplier = max(max_multiplier, pass_multiplier)
            max_abs_coefficient = max(
                max_abs_coefficient,
                pass_max_coefficient,
            )

            mean_estimate = statistics.fmean(estimates)
            if pass_index == 1:
                empirical_sd: float | str = ""
                empirical_stderr: float | str = ""
            else:
                empirical_sd = statistics.stdev(estimates)
                empirical_stderr = empirical_sd / math.sqrt(pass_index)

            progressive_rows.append(
                {
                    "family": run["family"],
                    "case": run["case"],
                    "tail_ratio": "1.00",
                    "passes": pass_index,
                    "seed": run["seed"],
                    "pass_estimate": f"{estimate:.17g}",
                    "pass_absolute_error": f"{abs(estimate - reference):.17g}",
                    "mean_estimate": f"{mean_estimate:.17g}",
                    "mean_absolute_error": f"{abs(mean_estimate - reference):.17g}",
                    "empirical_sd": (
                        "" if empirical_sd == "" else f"{empirical_sd:.17g}"
                    ),
                    "empirical_stderr": (
                        ""
                        if empirical_stderr == ""
                        else f"{empirical_stderr:.17g}"
                    ),
                    "pass_runtime_s": f"{runtime:.9f}",
                    "mean_runtime_s": f"{statistics.fmean(runtimes):.9f}",
                    "total_runtime_s": f"{sum(runtimes):.9f}",
                    "max_peak_rss_mb": f"{max(peak_rss_values):.3f}",
                    "pass_max_importance_multiplier": f"{pass_multiplier:.17g}",
                    "max_importance_multiplier": f"{max_multiplier:.17g}",
                    "pass_max_post_truncation_abs_coefficient": (
                        f"{pass_max_coefficient:.17g}"
                    ),
                    "max_post_truncation_abs_coefficient": (
                        f"{max_abs_coefficient:.17g}"
                    ),
                    "reference": f"{reference:.17g}",
                }
            )

            if pass_index in {10, 20, 50, max_passes}:
                sd_text = (
                    "n/a" if empirical_sd == "" else f"{empirical_sd:.3e}"
                )
                print(
                    f"case {case_index}, passes {pass_index}: "
                    f"error {abs(mean_estimate - reference):.3e}, "
                    f"SD {sd_text}, max weight {max_multiplier:.2e}",
                    flush=True,
                )

    single_path.parent.mkdir(parents=True, exist_ok=True)
    progressive_path.parent.mkdir(parents=True, exist_ok=True)
    with single_path.open("w", newline="") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=SINGLE_PASS_COLUMNS)
        writer.writeheader()
        writer.writerows(single_rows)
    with progressive_path.open("w", newline="") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=PROGRESSIVE_COLUMNS)
        writer.writeheader()
        writer.writerows(progressive_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
