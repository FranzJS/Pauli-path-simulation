#include "pauli_bench/truncation.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace pauli_bench {

void validate_k_strategy(
    const KStrategyConfig& config,
    bool allow_l1_mass) {
    switch (config.strategy) {
        case KStrategy::L1Mass:
            if (!allow_l1_mass) {
                throw std::invalid_argument(
                    "L1 K determination is only valid for deterministic propagation");
            }
            if (!std::isfinite(config.l1_cutoff) ||
                config.l1_cutoff < 0.0 || config.l1_cutoff > 1.0) {
                throw std::invalid_argument(
                    "L1 cutoff must be finite and in [0, 1]");
            }
            break;
        case KStrategy::SupportBudget:
            if (!std::isfinite(config.minimum_magnitude) ||
                config.minimum_magnitude < 0.0) {
                throw std::invalid_argument(
                    "minimum coefficient magnitude must be finite and nonnegative");
            }
            break;
        case KStrategy::Schedule:
            break;
    }
}

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

KDecision determine_k(
    Frontier& frontier,
    const KStrategyConfig& config,
    std::size_t truncation_index) {
    validate_k_strategy(config);
    switch (config.strategy) {
        case KStrategy::L1Mass: {
            const auto partition = partition_l1_mass(
                frontier, config.l1_cutoff);
            return {
                frontier.size() - partition.tail_terms,
                partition.total_l1_mass > 0.0 && config.l1_cutoff > 0.0};
        }
        case KStrategy::SupportBudget:
            return {
                support_budget_size(
                    frontier,
                    config.maximum_support,
                    config.minimum_magnitude),
                false};
        case KStrategy::Schedule:
            if (truncation_index >= config.schedule.size()) {
                throw std::invalid_argument("K schedule has too few events");
            }
            return {config.schedule[truncation_index], false};
    }
    throw std::invalid_argument("unknown K determination strategy");
}

std::size_t retain_largest_k(
    Frontier& frontier,
    std::size_t retained_terms,
    bool already_sorted_by_magnitude) {
    retained_terms = std::min(retained_terms, frontier.size());
    const std::size_t removed = frontier.size() - retained_terms;
    if (removed == 0) {
        return 0;
    }
    if (!already_sorted_by_magnitude) {
        sort_frontier_by_magnitude(frontier);
    }
    std::move(
        frontier.begin() + static_cast<std::ptrdiff_t>(removed),
        frontier.end(),
        frontier.begin());
    frontier.resize(retained_terms);
    return removed;
}

std::size_t truncate_l1_mass(
    Frontier& frontier,
    double cutoff_fraction) {
    KStrategyConfig config;
    config.strategy = KStrategy::L1Mass;
    config.l1_cutoff = cutoff_fraction;
    const auto decision = determine_k(frontier, config, 0);
    return retain_largest_k(
        frontier,
        decision.retained_terms,
        decision.frontier_sorted_by_magnitude);
}

std::size_t support_budget_size(
    const Frontier& frontier,
    std::size_t maximum_support,
    double minimum_magnitude) {
    KStrategyConfig config;
    config.strategy = KStrategy::SupportBudget;
    config.maximum_support = maximum_support;
    config.minimum_magnitude = minimum_magnitude;
    validate_k_strategy(config);
    const auto above_threshold = static_cast<std::size_t>(std::count_if(
        frontier.begin(),
        frontier.end(),
        [&](const Term& term) {
            return std::abs(term.coefficient) >= minimum_magnitude;
        }));
    return std::min(maximum_support, above_threshold);
}

std::size_t truncate_support_budget(
    Frontier& frontier,
    std::size_t maximum_support,
    double minimum_magnitude) {
    KStrategyConfig config;
    config.strategy = KStrategy::SupportBudget;
    config.maximum_support = maximum_support;
    config.minimum_magnitude = minimum_magnitude;
    const auto decision = determine_k(frontier, config, 0);
    return retain_largest_k(frontier, decision.retained_terms);
}

}  // namespace pauli_bench
