#pragma once

#include "pauli_bench/types.hpp"

#include <cstdint>
#include <vector>

namespace pauli_bench {

Circuit make_clifford_t_brickwork(
    int qubits,
    int layers,
    double t_density,
    std::uint64_t seed,
    double depolarizing_probability = 0.0);

Circuit make_nonintegrable_ising(
    int qubits,
    int steps,
    double dt,
    double coupling,
    double transverse_field,
    double longitudinal_field);

std::vector<Circuit> benchmark_circuits();

}  // namespace pauli_bench
