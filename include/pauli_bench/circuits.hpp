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

// Build U followed by its exact gate-by-gate inverse. `forward_layers` is the
// depth of U, so the returned circuit has 2 * forward_layers logical layers.
// Supported base models are the noiseless "clifford_t" and "ising" models.
Circuit make_identity_echo_circuit(
    int qubits,
    int forward_layers,
    const std::string& base_model,
    double l1_cutoff = -1.0);

// Construct W followed by U and then the exact gate-by-gate inverse of U.
// A single Clifford seed reproducibly generates W and U; U retains the same
// realization as the corresponding unprefixed circuit with that seed.
Circuit make_prefixed_identity_echo_circuit(
    int qubits,
    int forward_layers,
    int prefix_layers,
    const std::string& base_model,
    std::uint64_t circuit_seed,
    double l1_cutoff = -1.0);

struct ExactIdentityReference {
    double value{};
    int lightcone_qubits{};
};

// Compute <0|W^dagger O W|0> on an exact reduced statevector containing only
// the backward lightcone of the centered target observable.
ExactIdentityReference compute_prefixed_identity_reference(
    int qubits,
    int prefix_layers,
    const std::string& base_model,
    std::uint64_t circuit_seed);

// Construct a scalable simulation using the benchmark's fixed model-specific
// physical parameters. In addition to the three standard models, the
// "clifford_t_identity" and "ising_identity" variants construct U followed
// by U^dagger. The observable is centered in the qubit register.
Circuit make_configured_circuit(
    int qubits,
    int layers,
    const std::string& model,
    double l1_cutoff);

std::vector<Circuit> benchmark_circuits();

}  // namespace pauli_bench
