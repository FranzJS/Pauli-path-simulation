#pragma once

#include "pauli_bench/types.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace pauli_bench {

Circuit make_clifford_t_brickwork(
    int qubits,
    int layers,
    double t_density,
    std::uint64_t seed,
    double depolarizing_probability = 0.0,
    double l1_cutoff = -1.0);

Circuit make_nonintegrable_ising(
    int qubits,
    int steps,
    double dt,
    double coupling,
    double transverse_field,
    double longitudinal_field,
    double l1_cutoff = -1.0);

// Construct a scalable simulation using the benchmark's fixed model-specific
// physical parameters. Supported models are "clifford_t", "ising", and
// "clifford_t_depol". The observable is centered in the qubit register.
Circuit make_configured_circuit(
    int qubits,
    int layers,
    const std::string& model,
    double l1_cutoff);

std::vector<Circuit> benchmark_circuits();

}  // namespace pauli_bench
