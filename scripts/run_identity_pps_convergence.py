#!/usr/bin/env python3
"""Plot PPS convergence for Clifford+T and Ising W/U/U-dagger circuits."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from run_pps_convergence import (
    default_workers,
    get_schedule,
    render_svg,
    run_result,
)
from validate_reference_registry import configuration_key, validate_registry


OUTPUT_COLUMNS = [
    "family", "case", "qubits", "w_layers", "forward_layers", "total_layers",
    "circuit_seed",
    "k_strategy", "l1_cutoff", "maximum_support", "minimum_magnitude",
    "batch", "passes", "seed", "pass_estimate", "running_mean_estimate",
    "running_mean_absolute_error",
    "average_batch_running_mean_absolute_error", "deterministic_estimate",
    "deterministic_absolute_error", "reference", "reference_method",
    "reference_id", "pass_runtime_s", "cumulative_pps_runtime_s", "schedule",
]

CASES = [
    ("clifford_t_identity", "clifford_t", "Clifford+T prefixed echo"),
    ("ising_identity", "ising", "Nonintegrable Ising prefixed echo"),
]

DEFAULT_REGISTRY = (
    Path(__file__).resolve().parents[1]
    / "references"
    / "reference_registry.json"
)


def seed_for(case_index: int, batch_index: int, pass_index: int) -> int:
    return (
        202708050000
        + 1_000_003 * batch_index
        + 1009 * case_index
        + pass_index
    )


def centered_observable(qubits: int, base_model: str) -> dict[str, str]:
    if base_model == "clifford_t":
        start = min(max(qubits // 2 - 3, 0), qubits - 4)
        x_mask = (1 << start) | (1 << (start + 2)) | (1 << (start + 3))
        z_mask = 1 << (start + 1)
    else:
        start = min(max(qubits // 2 - 1, 0), qubits - 2)
        x_mask = 0
        z_mask = (1 << start) | (1 << (start + 1))
    return {"x_mask": hex(x_mask), "z_mask": hex(z_mask)}


def reference_entry(
    qubits: int,
    forward_layers: int,
    w_layers: int,
    model: str,
    base_model: str,
    circuit_seed: int,
) -> dict[str, Any]:
    is_clifford = base_model == "clifford_t"
    identifier = (
        f"{model}_n{qubits}_u{forward_layers}_w{w_layers}_v2"
        + (f"_seed{circuit_seed}" if is_clifford else "")
    )
    return {
        "id": identifier,
        "model": model,
        "qubits": qubits,
        "layers": forward_layers,
        "circuit_generation_version": 2,
        "circuit_seed": circuit_seed if is_clifford else None,
        "parameters": {
            "t_density": 0.70 if is_clifford else None,
            "depolarizing_probability": 0.0 if is_clifford else None,
            "dt": None if is_clifford else 0.12,
            "coupling": None if is_clifford else 1.0,
            "transverse_field": None if is_clifford else 0.91,
            "longitudinal_field": None if is_clifford else 0.37,
            "prefix_depth": w_layers,
        },
        "observable": centered_observable(qubits, base_model),
        "reference": {
            "value": 0.0,
            "method": "lightcone_statevector",
            "precision": "float64",
            "uncertainty": None,
        },
        "provenance": {
            "software": "pauli_bench",
            "software_version": "circuit_generation_v2",
            "generated_at": None,
            "notes": "",
        },
    }


def ensure_reference(
    registry_path: Path,
    reference_executable: Path,
    qubits: int,
    forward_layers: int,
    w_layers: int,
    model: str,
    base_model: str,
    circuit_seed: int,
) -> dict[str, Any]:
    validate_registry(registry_path, quiet=True)
    with registry_path.open(encoding="utf-8") as input_file:
        registry = json.load(input_file)
    if registry.get("schema_version") != 2:
        raise RuntimeError(
            f"{registry_path} must use reference registry schema version 2"
        )

    candidate = reference_entry(
        qubits,
        forward_layers,
        w_layers,
        model,
        base_model,
        circuit_seed,
    )
    candidate_key = configuration_key(candidate)
    for entry in registry.get("references", []):
        if configuration_key(entry) == candidate_key:
            print(f"using cached exact reference {entry['id']}", flush=True)
            return entry

    completed = subprocess.run(
        [
            str(reference_executable),
            str(qubits),
            str(w_layers),
            base_model,
            str(circuit_seed),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    fields = next(csv.reader([completed.stdout.strip()]))
    if len(fields) != 2:
        raise RuntimeError(
            f"unexpected exact-reference output: {completed.stdout!r}"
        )
    value = float(fields[0])
    lightcone_qubits = int(fields[1])
    if not math.isfinite(value):
        raise RuntimeError("exact reference is not finite")

    candidate["reference"]["value"] = value
    candidate["provenance"]["generated_at"] = datetime.now(
        timezone.utc
    ).isoformat()
    candidate["provenance"]["notes"] = (
        f"Exact statevector on a {lightcone_qubits}-qubit backward lightcone "
        f"of the depth-{w_layers} W circuit."
    )
    registry["references"].append(candidate)

    registry_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=registry_path.parent,
            prefix=f".{registry_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as output_file:
            json.dump(registry, output_file, indent=2)
            output_file.write("\n")
            temporary_path = Path(output_file.name)
        validate_registry(temporary_path, quiet=True)
        os.chmod(temporary_path, registry_path.stat().st_mode & 0o777)
        os.replace(temporary_path, registry_path)
    except BaseException:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
        raise

    print(
        f"computed and registered exact reference {candidate['id']}: "
        f"{value:.17g} ({lightcone_qubits}-qubit lightcone)",
        flush=True,
    )
    return candidate


def parse_l1_cutoffs(encoded: str) -> tuple[float, float]:
    values = tuple(float(value) for value in encoded.split(","))
    if len(values) == 1:
        values = values * 2
    if len(values) != 2:
        raise ValueError("L1 cutoffs must contain one value or two values")
    if any(
        not math.isfinite(value) or value < 0.0 or value > 1.0
        for value in values
    ):
        raise ValueError("every L1 cutoff must be finite and in [0, 1]")
    return values


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run scheduled PPS convergence on n-qubit circuits made from a "
            "shallow W, n U layers, and the exact inverse of U."
        )
    )
    parser.add_argument("build_dir", type=Path)
    parser.add_argument("qubits", type=int)
    subparsers = parser.add_subparsers(dest="strategy", required=True)

    support = subparsers.add_parser("support")
    support.add_argument("maximum_support", type=int)
    support.add_argument("minimum_magnitude", type=float)

    l1 = subparsers.add_parser("l1")
    l1.add_argument(
        "cutoffs",
        nargs="?",
        default="0.00625,0.00005",
        help=(
            "one shared cutoff or Clifford+T,Ising cutoffs "
            "(default: 0.00625,0.00005)"
        ),
    )

    for subparser in (support, l1):
        subparser.add_argument(
            "--w-depth",
            type=int,
            default=5,
            help="depth of the initial W circuit (default: 5)",
        )
        subparser.add_argument(
            "--circuit-seed",
            type=int,
            default=20260715,
            help="master seed reproducibly selecting Clifford+T W and U",
        )
        subparser.add_argument(
            "--registry",
            type=Path,
            default=DEFAULT_REGISTRY,
            help="JSON reference registry to read and update",
        )
        subparser.add_argument("--passes", type=int, default=100)
        subparser.add_argument("--batches", type=int, default=1)
        subparser.add_argument("--output", type=Path)
        subparser.add_argument("--workers", type=int, default=default_workers())
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    build_dir: Path = arguments.build_dir
    qubits: int = arguments.qubits
    strategy: str = arguments.strategy
    passes: int = arguments.passes
    batches: int = arguments.batches
    workers: int = arguments.workers
    w_layers: int = arguments.w_depth
    circuit_seed: int = arguments.circuit_seed
    registry_path: Path = arguments.registry.resolve()

    if qubits < 4 or qubits > 64:
        raise ValueError(
            "qubits must be in [4, 64] so both circuit families are valid"
        )
    if passes <= 0:
        raise ValueError("passes must be positive")
    if batches <= 0:
        raise ValueError("batches must be positive")
    if workers <= 0:
        raise ValueError("workers must be positive")
    if w_layers < 0:
        raise ValueError("W depth must be nonnegative")
    if circuit_seed < 0:
        raise ValueError("circuit seed must be nonnegative")

    if strategy == "support":
        maximum_support: int = arguments.maximum_support
        minimum_magnitude: float = arguments.minimum_magnitude
        if maximum_support < 0:
            raise ValueError("maximum support must be nonnegative")
        if not math.isfinite(minimum_magnitude) or minimum_magnitude < 0.0:
            raise ValueError("minimum magnitude must be finite and nonnegative")
        l1_cutoffs = (math.nan, math.nan)
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

    default_suffix = f"{batches}x{passes}" if batches > 1 else str(passes)
    output_prefix: Path = arguments.output or Path(
        f"results/identity_{strategy}_n{qubits}_w{w_layers}_{default_suffix}"
    )
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    executable = build_dir / "pauli_magnitude_pps_benchmark"
    reference_executable = build_dir / "pauli_identity_reference"
    if not executable.is_file():
        raise FileNotFoundError(
            f"benchmark executable not found: {executable}; build it first"
        )
    if not reference_executable.is_file():
        raise FileNotFoundError(
            f"reference executable not found: {reference_executable}; "
            "build it first"
        )
    if not registry_path.is_file():
        raise FileNotFoundError(f"reference registry not found: {registry_path}")
    os.environ["PAULI_REFERENCE_REGISTRY"] = str(registry_path)

    forward_layers = qubits
    total_layers = w_layers + 2 * qubits
    rows: list[dict[str, object]] = []
    trajectories: list[dict[str, object]] = []

    for case_index, (model, base_model, title) in enumerate(CASES):
        stored_reference = ensure_reference(
            registry_path,
            reference_executable,
            qubits,
            forward_layers,
            w_layers,
            model,
            base_model,
            circuit_seed,
        )
        exact_reference = float(stored_reference["reference"]["value"])
        circuit_arguments = [
            "prefix",
            str(w_layers),
            str(circuit_seed),
        ]
        l1_cutoff = l1_cutoffs[case_index]
        k_arguments = (
            ["support", str(maximum_support), str(minimum_magnitude)]
            if strategy == "support"
            else ["l1", str(l1_cutoff)]
        )
        schedule = get_schedule(
            executable,
            qubits,
            forward_layers,
            model,
            k_arguments,
            circuit_arguments,
        )
        encoded_schedule = ":".join(str(value) for value in schedule)
        deterministic = run_result([
            str(executable), str(qubits), str(forward_layers), model,
            *circuit_arguments,
            *k_arguments, "bfs",
        ])
        reported_reference = float(deterministic["reference"])
        if reported_reference != exact_reference:
            raise RuntimeError(
                f"{model} reported reference {reported_reference}, expected "
                f"stored W reference {exact_reference}"
            )
        if deterministic["reference_method"] != "lightcone_statevector":
            raise RuntimeError(
                f"{model} did not load its lightcone statevector reference"
            )
        if deterministic["reference_id"] != stored_reference["id"]:
            raise RuntimeError(
                f"{model} loaded reference {deterministic['reference_id']}, "
                f"expected {stored_reference['id']}"
            )

        deterministic_estimate = float(deterministic["estimate"])
        deterministic_error = abs(deterministic_estimate - exact_reference)
        print(
            f"{title}, n={qubits}, BFS: absolute error "
            f"{deterministic_error:.3e}",
            flush=True,
        )

        def run_pass(task: tuple[int, int]) -> dict[str, str]:
            batch_index, pass_index = task
            return run_result([
                str(executable), str(qubits), str(forward_layers), model,
                *circuit_arguments,
                "schedule", encoded_schedule, "pps_ht",
                str(seed_for(case_index, batch_index, pass_index)),
            ])

        tasks = [
            (batch_index, pass_index)
            for batch_index in range(batches)
            for pass_index in range(1, passes + 1)
        ]
        with ThreadPoolExecutor(max_workers=min(workers, len(tasks))) as executor:
            ordered_runs = list(executor.map(run_pass, tasks))

        batch_records: list[list[dict[str, object]]] = []
        batch_errors: list[list[float]] = []
        for batch_index in range(batches):
            batch_runs = ordered_runs[
                batch_index * passes:(batch_index + 1) * passes
            ]
            estimates: list[float] = []
            running_errors: list[float] = []
            cumulative_runtime = 0.0
            records: list[dict[str, object]] = []
            for pass_index, run in enumerate(batch_runs, start=1):
                seed = seed_for(case_index, batch_index, pass_index)
                estimate = float(run["estimate"])
                runtime = float(run["runtime_s"])
                estimates.append(estimate)
                cumulative_runtime += runtime
                running_mean = statistics.fmean(estimates)
                running_error = abs(running_mean - exact_reference)
                running_errors.append(running_error)
                records.append({
                    "family": run["family"],
                    "case": run["case"],
                    "qubits": qubits,
                    "w_layers": w_layers,
                    "forward_layers": forward_layers,
                    "total_layers": total_layers,
                    "circuit_seed": (
                        circuit_seed if base_model == "clifford_t" else ""
                    ),
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
                    "batch": batch_index + 1,
                    "passes": pass_index,
                    "seed": seed,
                    "pass_estimate": f"{estimate:.17g}",
                    "running_mean_estimate": f"{running_mean:.17g}",
                    "running_mean_absolute_error": f"{running_error:.17g}",
                    "deterministic_estimate": f"{deterministic_estimate:.17g}",
                    "deterministic_absolute_error": f"{deterministic_error:.17g}",
                    "reference": f"{exact_reference:.17g}",
                    "reference_method": deterministic["reference_method"],
                    "reference_id": deterministic["reference_id"],
                    "pass_runtime_s": f"{runtime:.9f}",
                    "cumulative_pps_runtime_s": f"{cumulative_runtime:.9f}",
                    "schedule": encoded_schedule,
                })
            batch_records.append(records)
            batch_errors.append(running_errors)
            print(
                f"{title}, batch {batch_index + 1}/{batches}, R={passes}: "
                f"running-mean error {running_errors[-1]:.3e}",
                flush=True,
            )

        average_errors = [
            statistics.fmean(errors[pass_index] for errors in batch_errors)
            for pass_index in range(passes)
        ]
        for records in batch_records:
            for pass_index, record in enumerate(records):
                record["average_batch_running_mean_absolute_error"] = (
                    f"{average_errors[pass_index]:.17g}"
                )
                rows.append(record)
        for pass_index in sorted({1, 10, 20, 50, passes}):
            if pass_index <= passes:
                print(
                    f"{title}, batch-average R={pass_index}: error "
                    f"{average_errors[pass_index - 1]:.3e}",
                    flush=True,
                )

        trajectories.append({
            "title": title,
            "batch_errors": batch_errors,
            "average_errors": average_errors,
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
        batches,
        strategy_description,
        deterministic_label,
        caption_lines=(
            f"n={qubits}; depth-{w_layers} W, n U layers, and n inverse-U "
            f"layers; {total_layers} total logical layers.",
            f"Clifford+T: density 0.70, circuit seed {circuit_seed}. Ising: "
            "dt=0.12, J=1, hx=0.91, hz=0.37. References are exact W "
            "lightcone statevectors.",
        ),
    )
    print(f"wrote {csv_path} and {output_prefix.with_suffix('.svg')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
