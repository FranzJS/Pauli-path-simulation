#pragma once

#include "pauli_bench/types.hpp"

#include <cstddef>

namespace pauli_bench {

struct L1MassPartition {
    std::size_t tail_terms{};
    double total_l1_mass{};
    double tail_l1_mass{};
};

void sort_frontier_by_magnitude(Frontier& frontier);

// Sorts the frontier by increasing magnitude, with the Pauli key as a
// deterministic tie-breaker, and identifies the largest removable prefix
// whose L1 mass does not exceed cutoff_fraction * total_l1_mass.
L1MassPartition partition_l1_mass(Frontier& frontier, double cutoff_fraction);

std::size_t truncate_l1_mass(Frontier& frontier, double cutoff_fraction);

}  // namespace pauli_bench
