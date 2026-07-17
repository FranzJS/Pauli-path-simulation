#include "pauli_bench/truncation.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace pauli_bench {

void sort_frontier_by_magnitude(Frontier& frontier) {
    std::sort(
        frontier.begin(),
        frontier.end(),
        [](const Term& lhs, const Term& rhs) {
            const double lhs_magnitude = std::abs(lhs.coefficient);
            const double rhs_magnitude = std::abs(rhs.coefficient);
            if (lhs_magnitude != rhs_magnitude) {
                return lhs_magnitude < rhs_magnitude;
            }
            if (lhs.pauli.x != rhs.pauli.x) {
                return lhs.pauli.x < rhs.pauli.x;
            }
            return lhs.pauli.z < rhs.pauli.z;
        });
}

L1MassPartition partition_l1_mass(
    Frontier& frontier,
    double cutoff_fraction) {
    L1MassPartition partition;
    if (frontier.empty() || cutoff_fraction <= 0.0) {
        return partition;
    }

    for (const auto& term : frontier) {
        partition.total_l1_mass += std::abs(term.coefficient);
    }

    const double removal_budget = cutoff_fraction * partition.total_l1_mass;
    if (removal_budget <= 0.0) {
        return partition;
    }

    sort_frontier_by_magnitude(frontier);

    while (partition.tail_terms < frontier.size()) {
        const double next_mass =
            std::abs(frontier[partition.tail_terms].coefficient);
        if (partition.tail_l1_mass + next_mass > removal_budget) {
            break;
        }
        partition.tail_l1_mass += next_mass;
        ++partition.tail_terms;
    }
    return partition;
}

std::size_t truncate_l1_mass(
    Frontier& frontier,
    double cutoff_fraction) {
    const auto partition = partition_l1_mass(frontier, cutoff_fraction);
    if (partition.tail_terms == 0) {
        return 0;
    }

    std::move(
        frontier.begin() + static_cast<std::ptrdiff_t>(partition.tail_terms),
        frontier.end(),
        frontier.begin());
    frontier.resize(frontier.size() - partition.tail_terms);
    return partition.tail_terms;
}

}  // namespace pauli_bench
