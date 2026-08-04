#!/usr/bin/env python3
"""Plot running-mean PPS error using a deterministic BFS-derived K schedule."""

from __future__ import annotations

import argparse
import csv
import math
import os
import statistics
import subprocess
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

OUTPUT_COLUMNS = [
    "family", "case", "qubits", "layers", "k_strategy", "l1_cutoff",
    "maximum_support", "minimum_magnitude", "passes", "seed", "pass_estimate",
    "running_mean_estimate", "running_mean_absolute_error",
    "empirical_standard_error",
    "deterministic_estimate", "deterministic_absolute_error", "reference",
    "pass_runtime_s", "cumulative_pps_runtime_s", "schedule",
]

CASES = [
    (20, 12, "clifford_t", "Clifford+T"),
    (20, 12, "ising", "Nonintegrable Ising"),
    (20, 12, "clifford_t_depol", "Noisy Clifford+T"),
]


def parse_result(stdout: str) -> dict[str, str]:
    values = next(csv.reader([stdout.strip()]))
    if len(values) != len(RAW_COLUMNS):
        raise RuntimeError(f"unexpected benchmark output: {stdout!r}")
    return dict(zip(RAW_COLUMNS, values, strict=True))


def run_result(command: list[str]) -> dict[str, str]:
    completed = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
    )
    return parse_result(completed.stdout)


def get_schedule(
    executable: Path,
    qubits: int,
    layers: int,
    model: str,
    k_arguments: list[str],
) -> tuple[int, ...]:
    completed = subprocess.run(
        [
            str(executable), str(qubits), str(layers), model,
            *k_arguments, "schedule",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    encoded = completed.stdout.strip()
    if not encoded:
        raise RuntimeError("deterministic run produced an empty K schedule")
    return tuple(int(value) for value in encoded.split(":"))


def seed_for(case_index: int, pass_index: int) -> int:
    return 202707270000 + 1009 * case_index + pass_index


def default_workers() -> int:
    return os.cpu_count() or 1


def parse_l1_cutoffs(encoded: str) -> tuple[float, float, float]:
    values = tuple(float(value) for value in encoded.split(","))
    if len(values) == 1:
        values = values * 3
    if len(values) != 3:
        raise ValueError("L1 cutoffs must contain one value or three values")
    if any(
        not math.isfinite(value) or value < 0.0 or value > 1.0
        for value in values
    ):
        raise ValueError("every L1 cutoff must be finite and in [0, 1]")
    return values


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run scheduled magnitude-PPS convergence using a deterministic "
            "support- or L1-derived K schedule."
        )
    )
    parser.add_argument("build_dir", type=Path)
    subparsers = parser.add_subparsers(dest="strategy", required=True)

    support = subparsers.add_parser("support")
    support.add_argument("maximum_support", type=int)
    support.add_argument("minimum_magnitude", type=float)

    l1 = subparsers.add_parser("l1")
    l1.add_argument(
        "cutoffs",
        nargs="?",
        default="0.00625,0.00005,0.034",
        help="one shared cutoff or three comma-separated per-case cutoffs",
    )

    for subparser in (support, l1):
        subparser.add_argument("--passes", type=int, default=100)
        subparser.add_argument("--output", type=Path)
        subparser.add_argument("--workers", type=int, default=default_workers())
    return parser.parse_args()


def pass_axis_ticks(passes: int, target_intervals: int = 5) -> list[int]:
    if passes <= 1:
        return [1]
    raw_step = (passes - 1) / target_intervals
    magnitude = 10 ** math.floor(math.log10(raw_step))
    normalized = raw_step / magnitude
    multiplier = next(value for value in (1, 2, 5, 10) if value >= normalized)
    step = max(1, int(multiplier * magnitude))
    ticks = [1, *range(step, passes + 1, step)]
    if ticks[-1] != passes:
        ticks.append(passes)
    return sorted(set(ticks))


def render_svg(
    output_path: Path,
    trajectories: list[dict[str, object]],
    passes: int,
    strategy_description: str,
    deterministic_label: str,
) -> None:
    width = 1200
    height = 475
    panel_width = 370
    panel_gap = 22
    left = 62
    top = 42
    plot_width = 295
    plot_height = 285
    blue = "#1f77b4"
    orange = "#d97706"
    standard_error_start = max(2, math.ceil(math.sqrt(passes)))
    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#20242b}'
        '.tick{font-size:11px}.title{font-size:15px;font-weight:600}'
        '.label{font-size:13px}.legend{font-size:11px}</style>',
    ]
    panel_labels = ["(a)", "(b)", "(c)"]
    for panel_index, (panel, trajectory) in enumerate(
        zip(panel_labels, trajectories, strict=True)
    ):
        origin_x = panel_index * (panel_width + panel_gap)
        x0 = origin_x + left
        y0 = top
        errors = [
            max(float(value), 1e-18)
            for value in trajectory["running_errors"]  # type: ignore[index]
        ]
        standard_errors = [
            None if value is None else max(float(value), 1e-18)
            for value in trajectory["standard_errors"]  # type: ignore[index]
        ]
        baseline = max(float(trajectory["deterministic_error"]), 1e-18)
        visible_standard_errors = [
            value
            for index, value in enumerate(standard_errors, start=1)
            if index >= standard_error_start and value is not None
        ]
        all_values = [*errors, *visible_standard_errors, baseline]
        low_decade = math.floor(math.log10(min(all_values)))
        high_decade = math.ceil(math.log10(max(all_values)))
        if high_decade <= low_decade:
            high_decade = low_decade + 1

        def x_position(pass_index: int) -> float:
            if passes == 1:
                return x0 + plot_width / 2
            return x0 + (pass_index - 1) * plot_width / (passes - 1)

        def y_position(value: float) -> float:
            fraction = (
                math.log10(max(value, 1e-18)) - low_decade
            ) / (high_decade - low_decade)
            return y0 + plot_height * (1.0 - fraction)

        svg.append(
            f'<text class="title" x="{x0 + plot_width / 2:.1f}" y="22" '
            f'text-anchor="middle">{panel} {trajectory["title"]}</text>'
        )
        svg.append(
            f'<rect x="{x0}" y="{y0}" width="{plot_width}" '
            f'height="{plot_height}" fill="none" stroke="#30343b"/>'
        )
        for decade in range(low_decade, high_decade + 1):
            y = y_position(10.0**decade)
            svg.append(
                f'<line x1="{x0}" y1="{y:.2f}" x2="{x0 + plot_width}" '
                f'y2="{y:.2f}" stroke="#d8dee8" stroke-width="1"/>'
            )
            svg.append(
                f'<text class="tick" x="{x0 - 7}" y="{y + 4:.2f}" '
                f'text-anchor="end">10<tspan baseline-shift="super" '
                f'font-size="8">{decade}</tspan></text>'
            )
        for tick in pass_axis_ticks(passes):
            x = x_position(tick)
            svg.append(
                f'<line x1="{x:.2f}" y1="{y0}" x2="{x:.2f}" '
                f'y2="{y0 + plot_height}" stroke="#e4e8ef" stroke-width="1"/>'
            )
            svg.append(
                f'<text class="tick" x="{x:.2f}" y="{y0 + plot_height + 18}" '
                f'text-anchor="middle">{tick}</text>'
            )
        points = " ".join(
            f"{x_position(index):.2f},{y_position(value):.2f}"
            for index, value in enumerate(errors, start=1)
        )
        svg.append(
            f'<polyline points="{points}" fill="none" stroke="{blue}" '
            f'stroke-width="1.6" stroke-linejoin="round"/>'
        )
        standard_error_points = " ".join(
            f"{x_position(index):.2f},{y_position(value):.2f}"
            for index, value in enumerate(standard_errors, start=1)
            if index >= standard_error_start and value is not None
        )
        if standard_error_points:
            svg.append(
                f'<polyline points="{standard_error_points}" fill="none" '
                f'stroke="{orange}" stroke-width="1.5" '
                f'stroke-dasharray="2,4" stroke-linejoin="round"/>'
            )
        baseline_y = y_position(baseline)
        svg.append(
            f'<line x1="{x0}" y1="{baseline_y:.2f}" '
            f'x2="{x0 + plot_width}" y2="{baseline_y:.2f}" '
            f'stroke="{blue}" stroke-width="1.4" stroke-dasharray="7,5"/>'
        )
        legend_x = x0 + plot_width - 119
        legend_y = y0 + 16
        svg.append(
            f'<rect x="{legend_x - 7}" y="{legend_y - 13}" width="126" '
            f'height="57" fill="white" fill-opacity="0.86"/>'
        )
        svg.append(
            f'<line x1="{legend_x}" y1="{legend_y}" x2="{legend_x + 30}" '
            f'y2="{legend_y}" stroke="{blue}" stroke-width="1.6"/>'
        )
        svg.append(
            f'<text class="legend" x="{legend_x + 35}" y="{legend_y + 4}">'
            'Magnitude PPS</text>'
        )
        svg.append(
            f'<line x1="{legend_x}" y1="{legend_y + 18}" '
            f'x2="{legend_x + 30}" y2="{legend_y + 18}" stroke="{blue}" '
            f'stroke-width="1.4" stroke-dasharray="7,5"/>'
        )
        svg.append(
            f'<text class="legend" x="{legend_x + 35}" y="{legend_y + 22}">'
            f'{deterministic_label}</text>'
        )
        svg.append(
            f'<line x1="{legend_x}" y1="{legend_y + 36}" '
            f'x2="{legend_x + 30}" y2="{legend_y + 36}" stroke="{orange}" '
            f'stroke-width="1.5" stroke-dasharray="2,4"/>'
        )
        svg.append(
            f'<text class="legend" x="{legend_x + 35}" y="{legend_y + 40}">'
            'PPS empirical SE</text>'
        )
        svg.append(
            f'<text class="label" x="{x0 + plot_width / 2:.1f}" '
            f'y="{y0 + plot_height + 44}" text-anchor="middle">'
            'Independent passes R</text>'
        )
        if panel_index == 0:
            center_y = y0 + plot_height / 2
            svg.append(
                f'<text class="label" x="16" y="{center_y:.1f}" '
                f'text-anchor="middle" transform="rotate(-90 16 {center_y:.1f})">'
                'Running-mean absolute error</text>'
            )
    svg.append(
        f'<text class="label" x="600" y="405" text-anchor="middle">'
        f'{passes} independent magnitude-PPS/HT passes; K schedule from '
        f'deterministic {strategy_description}; truncation every 4 RZ gates.'
        '</text>'
    )
    svg.append(
        '<text class="tick" x="600" y="428" text-anchor="middle">'
        'All cases: 20 qubits and 12 layers/steps. Clifford+T: T density 0.70, '
        'circuit seed 20260715; noisy case: depolarizing p=0.05 after each layer.'
        '</text>'
    )
    svg.append(
        f'<text class="tick" x="600" y="449" text-anchor="middle">'
        'Ising: dt=0.12, J=1.0, hx=0.91, hz=0.37. Solid: running-mean error; '
        f'dashed: {deterministic_label} error.'
        '</text>'
    )
    svg.append('</svg>')
    output_path.write_text("\n".join(svg) + "\n", encoding="utf-8")


def main() -> int:
    arguments = parse_arguments()
    build_dir: Path = arguments.build_dir
    strategy: str = arguments.strategy
    passes: int = arguments.passes
    workers: int = arguments.workers
    if strategy == "support":
        maximum_support: int = arguments.maximum_support
        minimum_magnitude: float = arguments.minimum_magnitude
        if maximum_support < 0:
            raise ValueError("maximum support must be nonnegative")
        if not math.isfinite(minimum_magnitude) or minimum_magnitude < 0.0:
            raise ValueError("minimum magnitude must be finite and nonnegative")
        l1_cutoffs = (math.nan, math.nan, math.nan)
        strategy_description = (
            f"support BFS; MAX_SUPPORT={maximum_support}; "
            f"MIN_MAGNITUDE={minimum_magnitude:.1e}"
        )
        deterministic_label = "BFS + support"
    else:
        maximum_support = 0
        minimum_magnitude = 0.0
        l1_cutoffs = parse_l1_cutoffs(arguments.cutoffs)
        strategy_description = (
            "L1 BFS; per-case cutoffs="
            + ",".join(f"{value:.5g}" for value in l1_cutoffs)
        )
        deterministic_label = "BFS + L1"
    output_prefix: Path = arguments.output or Path(
        f"results/{strategy}_pps_convergence_{passes}"
    )
    if passes <= 0:
        raise ValueError("passes must be positive")
    if workers <= 0:
        raise ValueError("workers must be positive")

    executable = build_dir / "pauli_magnitude_pps_benchmark"
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    trajectories: list[dict[str, object]] = []

    for case_index, (qubits, layers, model, title) in enumerate(CASES):
        l1_cutoff = l1_cutoffs[case_index]
        k_arguments = (
            ["support", str(maximum_support), str(minimum_magnitude)]
            if strategy == "support"
            else ["l1", str(l1_cutoff)]
        )
        schedule = get_schedule(
            executable,
            qubits,
            layers,
            model,
            k_arguments,
        )
        encoded_schedule = ":".join(str(value) for value in schedule)
        deterministic = run_result([
            str(executable), str(qubits), str(layers), model,
            *k_arguments, "bfs",
        ])
        reference = float(deterministic["reference"])
        if not math.isfinite(reference):
            raise RuntimeError(f"case {model} has no finite stored reference")
        deterministic_estimate = float(deterministic["estimate"])
        deterministic_error = abs(deterministic_estimate - reference)
        print(
            f"{title}, BFS: absolute error {deterministic_error:.3e}",
            flush=True,
        )

        def run_pass(pass_index: int) -> dict[str, str]:
            return run_result([
                str(executable), str(qubits), str(layers), model,
                "schedule", encoded_schedule, "pps_ht",
                str(seed_for(case_index, pass_index)),
            ])

        worker_count = min(workers, passes)
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            ordered_runs = list(executor.map(run_pass, range(1, passes + 1)))

        estimates: list[float] = []
        running_errors: list[float] = []
        standard_errors: list[float | None] = []
        cumulative_runtime = 0.0
        for pass_index, run in enumerate(ordered_runs, start=1):
            seed = seed_for(case_index, pass_index)
            estimate = float(run["estimate"])
            runtime = float(run["runtime_s"])
            estimates.append(estimate)
            cumulative_runtime += runtime
            running_mean = sum(estimates) / pass_index
            running_error = abs(running_mean - reference)
            running_errors.append(running_error)
            standard_errors.append(
                None
                if pass_index < 2
                else statistics.stdev(estimates) / math.sqrt(pass_index)
            )
            standard_error = standard_errors[-1]
            rows.append({
                "family": run["family"],
                "case": run["case"],
                "qubits": qubits,
                "layers": layers,
                "k_strategy": strategy,
                "l1_cutoff": (
                    f"{l1_cutoff:.17g}" if strategy == "l1" else ""
                ),
                "maximum_support": (
                    maximum_support if strategy == "support" else ""
                ),
                "minimum_magnitude": (
                    f"{minimum_magnitude:.17g}"
                    if strategy == "support"
                    else ""
                ),
                "passes": pass_index,
                "seed": seed,
                "pass_estimate": f"{estimate:.17g}",
                "running_mean_estimate": f"{running_mean:.17g}",
                "running_mean_absolute_error": f"{running_error:.17g}",
                "empirical_standard_error": (
                    ""
                    if standard_error is None
                    else f"{standard_error:.17g}"
                ),
                "deterministic_estimate": f"{deterministic_estimate:.17g}",
                "deterministic_absolute_error": f"{deterministic_error:.17g}",
                "reference": f"{reference:.17g}",
                "pass_runtime_s": f"{runtime:.9f}",
                "cumulative_pps_runtime_s": f"{cumulative_runtime:.9f}",
                "schedule": encoded_schedule,
            })
            if pass_index in {1, 10, 20, 50, passes}:
                print(
                    f"{title}, R={pass_index}: running-mean error "
                    f"{running_error:.3e}",
                    flush=True,
                )

        trajectories.append({
            "title": title,
            "running_errors": running_errors,
            "standard_errors": standard_errors,
            "deterministic_error": deterministic_error,
        })

    csv_path = output_prefix.with_suffix(".csv")
    with csv_path.open("w", newline="") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    render_svg(
        output_prefix.with_suffix(".svg"),
        trajectories,
        passes,
        strategy_description,
        deterministic_label,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
