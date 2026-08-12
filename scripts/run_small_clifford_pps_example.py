#!/usr/bin/env python3
"""Reproduce a finite-sample PPS error reversal and coefficient staircase."""

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
    "reference_id",
    "absolute_error", "runtime_s", "peak_support_terms",
    "peak_pre_truncation_terms", "peak_post_truncation_terms",
    "truncation_events", "peak_vector_capacity_terms",
    "peak_merge_entries", "peak_merge_capacity_slots", "max_heavy_terms",
    "max_sampled_tail_terms", "total_sampled_tail_terms",
    "max_importance_multiplier", "max_post_truncation_abs_coefficient",
]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build_dir", type=Path)
    parser.add_argument("--qubits", type=int, default=10)
    parser.add_argument("--layers", type=int, default=12)
    parser.add_argument("--supports", default="10000,30000")
    parser.add_argument("--passes", type=int, default=80)
    parser.add_argument("--seed-base", type=int, default=81000000)
    parser.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/small_clifford_pps_counterexample"),
    )
    return parser.parse_args()


def run_result(command: list[str]) -> dict[str, str]:
    completed = subprocess.run(
        command, check=True, capture_output=True, text=True
    )
    values = next(csv.reader([completed.stdout.strip()]))
    if len(values) != len(RAW_COLUMNS):
        raise RuntimeError(f"unexpected benchmark output: {completed.stdout!r}")
    return dict(zip(RAW_COLUMNS, values, strict=True))


def load_spectrum(
    executable: Path, qubits: int, layers: int, minimum_terms: int
) -> tuple[float, int, int, list[tuple[int, int, float, int]]]:
    text = subprocess.run(
        [
            str(executable), str(qubits), str(layers), "clifford_t",
            str(minimum_terms),
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()
    reference = float(text[0].split(",", 1)[1])
    event = int(text[1].split(",", 1)[1])
    frontier_terms = int(text[2].split(",", 1)[1])
    plateaus = [
        (int(row[0]), int(row[1]), float(row[2]), int(row[3]))
        for row in csv.reader(text[4:])
    ]
    return reference, event, frontier_terms, plateaus


def render_svg(
    path: Path,
    trajectories: list[dict[str, object]],
    reference: float,
    event: int,
    frontier_terms: int,
    plateaus: list[tuple[int, int, float, int]],
    passes: int,
    qubits: int,
    layers: int,
) -> None:
    width, height = 1100, 490
    left_x, right_x, top = 75, 625, 48
    plot_w, plot_h = 400, 310
    colors = ["#1f77b4", "#d97706", "#2ca02c", "#9467bd"]
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#20242b}'
        '.tick{font-size:11px}.title{font-size:15px;font-weight:600}'
        '.label{font-size:13px}.legend{font-size:11px}</style>',
    ]

    all_errors = [
        max(value, 1e-18)
        for trajectory in trajectories
        for value in trajectory["running_errors"]  # type: ignore[index]
    ]
    all_errors += [
        max(value, 1e-18)
        for trajectory in trajectories
        for value in trajectory["standard_errors"]  # type: ignore[index]
        if value is not None
    ]
    y_low = math.floor(math.log10(min(all_errors)))
    y_high = math.ceil(math.log10(max(all_errors)))

    def error_x(index: int) -> float:
        return left_x + (index - 1) * plot_w / max(passes - 1, 1)

    def error_y(value: float) -> float:
        fraction = (math.log10(max(value, 1e-18)) - y_low) / (y_high - y_low)
        return top + plot_h * (1 - fraction)

    svg += [
        f'<text class="title" x="{left_x + plot_w / 2}" y="24" '
        'text-anchor="middle">(a) A larger K gives a worse realized mean</text>',
        f'<rect x="{left_x}" y="{top}" width="{plot_w}" height="{plot_h}" '
        'fill="none" stroke="#30343b"/>',
    ]
    for decade in range(y_low, y_high + 1):
        y = error_y(10.0**decade)
        svg.append(
            f'<line x1="{left_x}" y1="{y:.2f}" x2="{left_x + plot_w}" '
            f'y2="{y:.2f}" stroke="#d8dee8"/>'
        )
        svg.append(
            f'<text class="tick" x="{left_x - 8}" y="{y + 4:.2f}" '
            f'text-anchor="end">10^{decade}</text>'
        )
    for tick in (1, 20, 40, 60, passes):
        if tick > passes:
            continue
        x = error_x(tick)
        svg.append(
            f'<text class="tick" x="{x:.2f}" y="{top + plot_h + 18}" '
            f'text-anchor="middle">{tick}</text>'
        )
    for index, trajectory in enumerate(trajectories):
        color = colors[index]
        errors = trajectory["running_errors"]
        ses = trajectory["standard_errors"]
        error_points = " ".join(
            f"{error_x(i):.2f},{error_y(value):.2f}"
            for i, value in enumerate(errors, 1)  # type: ignore[arg-type]
        )
        se_points = " ".join(
            f"{error_x(i):.2f},{error_y(value):.2f}"
            for i, value in enumerate(ses, 1)  # type: ignore[arg-type]
            if value is not None and i >= 10
        )
        svg.append(
            f'<polyline points="{error_points}" fill="none" stroke="{color}" '
            'stroke-width="1.6"/>'
        )
        svg.append(
            f'<polyline points="{se_points}" fill="none" stroke="{color}" '
            'stroke-width="1.2" stroke-dasharray="3,4"/>'
        )
        legend_y = top + 17 + 37 * index
        svg.append(
            f'<line x1="{left_x + 238}" y1="{legend_y}" '
            f'x2="{left_x + 266}" y2="{legend_y}" stroke="{color}" '
            'stroke-width="1.6"/>'
        )
        svg.append(
            f'<text class="legend" x="{left_x + 272}" y="{legend_y + 4}">'
            f'K={trajectory["support"]:,}, error={trajectory["final_error"]:.2g}'
            '</text>'
        )
        svg.append(
            f'<line x1="{left_x + 238}" y1="{legend_y + 14}" '
            f'x2="{left_x + 266}" y2="{legend_y + 14}" stroke="{color}" '
            'stroke-width="1.2" stroke-dasharray="3,4"/>'
        )
    svg += [
        f'<text class="label" x="{left_x + plot_w / 2}" y="{top + plot_h + 44}" '
        'text-anchor="middle">Independent PPS passes R (dashed: empirical SE)</text>',
        f'<text class="label" x="18" y="{top + plot_h / 2}" '
        f'transform="rotate(-90 18 {top + plot_h / 2})" text-anchor="middle">'
        'Running-mean absolute error</text>',
    ]

    positive = [row for row in plateaus if row[2] > 0.0]
    mag_low = math.floor(math.log10(min(row[2] for row in positive)))
    mag_high = math.ceil(math.log10(max(row[2] for row in positive)))
    rank_high = math.log10(max(frontier_terms, 10))

    def spectrum_x(rank: int) -> float:
        return right_x + math.log10(max(rank, 1)) * plot_w / rank_high

    def spectrum_y(magnitude: float) -> float:
        fraction = (math.log10(magnitude) - mag_low) / (mag_high - mag_low)
        return top + plot_h * (1 - fraction)

    svg += [
        f'<text class="title" x="{right_x + plot_w / 2}" y="24" '
        'text-anchor="middle">(b) Exact pre-truncation coefficient staircase</text>',
        f'<rect x="{right_x}" y="{top}" width="{plot_w}" height="{plot_h}" '
        'fill="none" stroke="#30343b"/>',
    ]
    for decade in range(mag_low, mag_high + 1):
        y = spectrum_y(10.0**decade)
        svg.append(
            f'<line x1="{right_x}" y1="{y:.2f}" x2="{right_x + plot_w}" '
            f'y2="{y:.2f}" stroke="#d8dee8"/>'
        )
        svg.append(
            f'<text class="tick" x="{right_x - 8}" y="{y + 4:.2f}" '
            f'text-anchor="end">10^{decade}</text>'
        )
    for decade in range(0, math.ceil(rank_high) + 1):
        rank = 10**decade
        if rank > frontier_terms:
            continue
        x = spectrum_x(rank)
        svg.append(
            f'<text class="tick" x="{x:.2f}" y="{top + plot_h + 18}" '
            f'text-anchor="middle">10^{decade}</text>'
        )
    for index, trajectory in enumerate(trajectories):
        support = int(trajectory["support"])
        if support > frontier_terms:
            continue
        x = spectrum_x(support)
        svg.append(
            f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" '
            f'y2="{top + plot_h}" stroke="{colors[index]}" '
            'stroke-width="1.2" stroke-dasharray="5,4"/>'
        )
        svg.append(
            f'<text class="legend" x="{x - 4:.2f}" y="{top + 15}" '
            f'text-anchor="end" fill="{colors[index]}">K={support:,}</text>'
        )
    for start, end, magnitude, _ in positive:
        y = spectrum_y(magnitude)
        svg.append(
            f'<line x1="{spectrum_x(start):.2f}" y1="{y:.2f}" '
            f'x2="{spectrum_x(end):.2f}" y2="{y:.2f}" '
            'stroke="#30343b" stroke-width="1.5"/>'
        )
    svg += [
        f'<text class="label" x="{right_x + plot_w / 2}" '
        f'y="{top + plot_h + 44}" text-anchor="middle">'
        'Coefficient rank (descending magnitude)</text>',
        f'<text class="label" x="{right_x - 52}" y="{top + plot_h / 2}" '
        f'transform="rotate(-90 {right_x - 52} {top + plot_h / 2})" '
        'text-anchor="middle">|coefficient|</text>',
        f'<text class="tick" x="550" y="432" text-anchor="middle">'
        f'{qubits}-qubit, {layers}-layer fixed Clifford+T circuit; '
        f'exact reference={reference:.8g}; first exact frontier above both budgets '
        f'is event '
        f'{event + 1}: {frontier_terms:,} terms.</text>',
        '<text class="tick" x="550" y="454" text-anchor="middle">'
        'The reversal is a finite-sample outcome; the dashed SE and repeated '
        'independent batches are required to assess expected precision.</text>',
        '</svg>',
    ]
    path.write_text("\n".join(svg) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_arguments()
    supports = tuple(int(value) for value in args.supports.split(","))
    if len(supports) < 2 or any(value <= 0 for value in supports):
        raise ValueError("--supports requires at least two positive integers")
    if args.passes < 2 or args.workers <= 0:
        raise ValueError("passes must be at least two and workers positive")

    benchmark = args.build_dir / "pauli_magnitude_pps_benchmark"
    spectrum_executable = args.build_dir / "pauli_frontier_spectrum"
    reference, event, frontier_terms, plateaus = load_spectrum(
        spectrum_executable, args.qubits, args.layers, max(supports) + 1
    )
    trajectories: list[dict[str, object]] = []
    rows: list[dict[str, object]] = []

    for support in supports:
        schedule = subprocess.run(
            [
                str(benchmark), str(args.qubits), str(args.layers),
                "clifford_t", "support", str(support), "0", "schedule",
            ],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()

        def run_pass(index: int) -> dict[str, str]:
            return run_result([
                str(benchmark), str(args.qubits), str(args.layers),
                "clifford_t", "schedule", schedule, "pps_ht",
                str(args.seed_base + index),
            ])

        with ThreadPoolExecutor(
            max_workers=min(args.workers, args.passes)
        ) as executor:
            runs = list(executor.map(run_pass, range(args.passes)))

        estimates: list[float] = []
        running_errors: list[float] = []
        standard_errors: list[float | None] = []
        for index, run in enumerate(runs, 1):
            estimates.append(float(run["estimate"]))
            mean = statistics.mean(estimates)
            error = abs(mean - reference)
            se = None if index == 1 else statistics.stdev(estimates) / math.sqrt(index)
            running_errors.append(error)
            standard_errors.append(se)
            rows.append({
                "maximum_support": support,
                "pass": index,
                "seed": args.seed_base + index - 1,
                "estimate": f'{estimates[-1]:.17g}',
                "running_mean": f'{mean:.17g}',
                "running_error": f'{error:.17g}',
                "empirical_standard_error": "" if se is None else f'{se:.17g}',
                "reference": f'{reference:.17g}',
                "schedule": schedule,
            })
        trajectories.append({
            "support": support,
            "running_errors": running_errors,
            "standard_errors": standard_errors,
            "final_error": running_errors[-1],
        })
        print(
            f"K={support}: final error={running_errors[-1]:.6g}, "
            f"SE={standard_errors[-1]:.6g}",
            flush=True,
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.with_suffix(".csv").open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with args.output.with_name(args.output.name + "_spectrum.csv").open(
        "w", newline=""
    ) as output:
        writer = csv.writer(output)
        writer.writerow(["event", "frontier_terms", "rank_start", "rank_end", "magnitude", "count"])
        for start, end, magnitude, count in plateaus:
            writer.writerow([event, frontier_terms, start, end, f'{magnitude:.17g}', count])
    render_svg(
        args.output.with_suffix(".svg"), trajectories, reference, event,
        frontier_terms, plateaus, args.passes, args.qubits, args.layers
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
