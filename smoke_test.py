"""Minimal repository execution smoke test."""

from __future__ import annotations


def main() -> None:
    n_qubits = 10
    state_dimension = 1 << n_qubits
    print(f"GitHub execution smoke test: {n_qubits} qubits -> dimension {state_dimension}")


if __name__ == "__main__":
    main()
