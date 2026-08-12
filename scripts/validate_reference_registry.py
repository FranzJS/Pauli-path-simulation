#!/usr/bin/env python3
"""Validate the JSON reference registry without regenerating references."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any


PARAMETER_NAMES = {
    "t_density",
    "depolarizing_probability",
    "dt",
    "coupling",
    "transverse_field",
    "longitudinal_field",
    "prefix_depth",
}
REFERENCE_METHODS = {
    "statevector",
    "lightcone_statevector",
    "converged_pauli",
}
MODELS = {
    "clifford_t",
    "clifford_t_depol",
    "ising",
    "clifford_t_identity",
    "ising_identity",
}


def integer(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        require(key not in result, f"duplicate JSON object key: {key}")
        result[key] = value
    return result


def finite_number(value: Any, field: str) -> float:
    require(
        isinstance(value, (int, float)) and not isinstance(value, bool),
        f"{field} must be numeric",
    )
    result = float(value)
    require(math.isfinite(result), f"{field} must be finite")
    return result


def nullable_number(value: Any, field: str) -> float | None:
    return None if value is None else finite_number(value, field)


def configuration_key(entry: dict[str, Any]) -> str:
    configuration = {
        "model": entry["model"],
        "qubits": entry["qubits"],
        "layers": entry["layers"],
        "circuit_generation_version": entry["circuit_generation_version"],
        "circuit_seed": entry["circuit_seed"],
        "parameters": entry["parameters"],
        "observable": entry["observable"],
    }
    return json.dumps(configuration, sort_keys=True, separators=(",", ":"))


def validate_entry(entry: Any, index: int) -> tuple[str, str, float]:
    prefix = f"references[{index}]"
    require(isinstance(entry, dict), f"{prefix} must be an object")
    required = {
        "id", "model", "qubits", "layers", "circuit_generation_version",
        "circuit_seed", "parameters", "observable", "reference", "provenance",
    }
    require(set(entry) == required, f"{prefix} fields must be exactly {sorted(required)}")
    require(
        isinstance(entry["id"], str)
        and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._+-]*", entry["id"]) is not None,
        f"{prefix}.id must be a nonempty CSV-safe identifier",
    )
    require(entry["model"] in MODELS, f"{prefix}.model is unknown")
    require(
        integer(entry["qubits"]) and 1 <= entry["qubits"] <= 64,
        f"{prefix}.qubits must be in [1, 64]",
    )
    require(
        integer(entry["layers"]) and entry["layers"] > 0,
        f"{prefix}.layers must be positive",
    )
    require(
        integer(entry["circuit_generation_version"])
        and entry["circuit_generation_version"] > 0,
        f"{prefix}.circuit_generation_version must be positive",
    )
    require(
        entry["circuit_seed"] is None
        or (integer(entry["circuit_seed"]) and entry["circuit_seed"] >= 0),
        f"{prefix}.circuit_seed must be null or nonnegative",
    )

    parameters = entry["parameters"]
    require(isinstance(parameters, dict), f"{prefix}.parameters must be an object")
    require(
        set(parameters) == PARAMETER_NAMES,
        f"{prefix}.parameters must contain exactly {sorted(PARAMETER_NAMES)}",
    )
    parsed_parameters = {
        name: nullable_number(value, f"{prefix}.parameters.{name}")
        for name, value in parameters.items()
    }
    identity_model = entry["model"].endswith("_identity")
    clifford_model = entry["model"].startswith("clifford_t")
    if identity_model:
        require(
            entry["circuit_generation_version"] >= 2,
            f"{prefix} identity model needs circuit generation version >= 2",
        )
        require(
            integer(parameters["prefix_depth"])
            and parameters["prefix_depth"] >= 0,
            f"{prefix}.parameters.prefix_depth must be a nonnegative integer",
        )
    else:
        require(
            parsed_parameters["prefix_depth"] is None,
            f"{prefix} standard model must have null prefix_depth",
        )

    if clifford_model:
        require(entry["circuit_seed"] is not None, f"{prefix} needs a circuit seed")
        require(parsed_parameters["t_density"] is not None, f"{prefix} needs t_density")
        require(
            parsed_parameters["depolarizing_probability"] is not None,
            f"{prefix} needs depolarizing_probability",
        )
        require(
            all(
                parsed_parameters[name] is None
                for name in (
                    "dt", "coupling", "transverse_field", "longitudinal_field"
                )
            ),
            f"{prefix} has Ising-only parameters",
        )
        require(
            0.0 <= parsed_parameters["t_density"] <= 1.0,
            f"{prefix}.parameters.t_density must be in [0, 1]",
        )
        require(
            0.0 <= parsed_parameters["depolarizing_probability"] <= 1.0,
            f"{prefix}.parameters.depolarizing_probability must be in [0, 1]",
        )
        if entry["model"] == "clifford_t_depol":
            require(
                parsed_parameters["depolarizing_probability"] > 0.0,
                f"{prefix} noisy model must have positive depolarizing probability",
            )
        else:
            require(
                parsed_parameters["depolarizing_probability"] == 0.0,
                f"{prefix} noiseless model must have zero depolarizing probability",
            )
    else:
        require(entry["circuit_seed"] is None, f"{prefix} Ising seed must be null")
        require(
            parsed_parameters["t_density"] is None
            and parsed_parameters["depolarizing_probability"] is None,
            f"{prefix} has Clifford-only parameters",
        )
        for name in ("dt", "coupling", "transverse_field", "longitudinal_field"):
            require(parsed_parameters[name] is not None, f"{prefix} needs {name}")
        require(parsed_parameters["dt"] > 0.0, f"{prefix}.parameters.dt must be positive")

    observable = entry["observable"]
    require(isinstance(observable, dict), f"{prefix}.observable must be an object")
    require(set(observable) == {"x_mask", "z_mask"}, f"{prefix}.observable fields invalid")
    try:
        x_mask = int(observable["x_mask"], 0)
        z_mask = int(observable["z_mask"], 0)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{prefix}.observable masks must be integer strings") from error
    require(x_mask >= 0 and z_mask >= 0, f"{prefix}.observable masks must be nonnegative")
    used_mask = x_mask | z_mask
    require(
        used_mask == 0 or used_mask.bit_length() <= entry["qubits"],
        f"{prefix}.observable addresses a qubit outside the register",
    )

    reference = entry["reference"]
    require(isinstance(reference, dict), f"{prefix}.reference must be an object")
    require(
        set(reference) == {"value", "method", "precision", "uncertainty"},
        f"{prefix}.reference fields are invalid",
    )
    value = finite_number(reference.get("value"), f"{prefix}.reference.value")
    require(reference.get("method") in REFERENCE_METHODS, f"{prefix}.reference.method invalid")
    if identity_model:
        require(
            reference["method"] == "lightcone_statevector",
            f"{prefix} identity reference must use lightcone_statevector",
        )
    require(
        isinstance(reference.get("precision"), str) and reference["precision"],
        f"{prefix}.reference.precision is invalid",
    )
    uncertainty = nullable_number(
        reference.get("uncertainty"), f"{prefix}.reference.uncertainty"
    )
    require(uncertainty is None or uncertainty >= 0.0, f"{prefix} uncertainty is negative")

    provenance = entry["provenance"]
    require(isinstance(provenance, dict), f"{prefix}.provenance must be an object")
    provenance_fields = {"software", "software_version", "generated_at", "notes"}
    require(
        set(provenance) == provenance_fields,
        f"{prefix}.provenance fields are invalid",
    )
    for field in ("software", "software_version", "notes"):
        require(
            isinstance(provenance[field], str),
            f"{prefix}.provenance.{field} must be a string",
        )
    require(
        provenance["generated_at"] is None
        or isinstance(provenance["generated_at"], str),
        f"{prefix}.provenance.generated_at must be null or a string",
    )
    return entry["id"], configuration_key(entry), value


def validate_registry(path: Path, *, quiet: bool = False) -> int:
    with path.open(encoding="utf-8") as input_file:
        document = json.load(input_file, object_pairs_hook=unique_object)
    require(isinstance(document, dict), "registry root must be an object")
    require(
        set(document) == {"schema_version", "references"},
        "registry root fields must be schema_version and references",
    )
    require(
        integer(document.get("schema_version"))
        and document["schema_version"] == 2,
        "unsupported schema_version",
    )
    references = document.get("references")
    require(isinstance(references, list), "references must be an array")

    ids: set[str] = set()
    configurations: dict[str, tuple[str, float]] = {}
    for index, entry in enumerate(references):
        identifier, key, value = validate_entry(entry, index)
        require(identifier not in ids, f"duplicate reference id: {identifier}")
        ids.add(identifier)
        if key in configurations:
            prior_id, prior_value = configurations[key]
            require(
                prior_value == value,
                f"conflicting references for one configuration: {prior_id}, {identifier}",
            )
            raise ValueError(
                f"duplicate configuration under ids {prior_id} and {identifier}"
            )
        configurations[key] = (identifier, value)
    if not quiet:
        print(f"validated {len(references)} reference records in {path}")
    return len(references)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "registry",
        nargs="?",
        type=Path,
        default=Path("references/reference_registry.json"),
    )
    arguments = parser.parse_args()
    validate_registry(arguments.registry)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
