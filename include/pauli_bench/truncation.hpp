#pragma once

#include "pauli_bench/types.hpp"

#include <cstddef>
#include <vector>

namespace pauli_bench {

struct L1MassPartition {
    std::size_t tail_terms{};
    double total_l1_mass{};
    double tail_l1_mass{};
};

enum class KStrategy { L1Mass, SupportBudget, Schedule };

struct KStrategyConfig {
    KStrategy strategy{KStrategy::L1Mass};
    double l1_cutoff{};
    std::size_t maximum_support{};
    double minimum_magnitude{};
    std::vector<std::size_t> schedule;
};

struct KDecision {
    std::size_t retained_terms{};
    bool frontier_sorted_by_magnitude{};
};

void validate_k_strategy(
    const KStrategyConfig& config,
    bool allow_l1_mass = true);

// Determine K for one truncation event. L1 determination sorts the frontier as
// part of computing its removable prefix; the decision records this so the
// deterministic selector does not sort it a second time.
KDecision determine_k(
    Frontier& frontier,
    const KStrategyConfig& config,
    std::size_t truncation_index);

// Deterministically retain the K largest-magnitude terms. The optional flag is
// supplied from KDecision when L1 determination already sorted the frontier.
std::size_t retain_largest_k(
    Frontier& frontier,
    std::size_t retained_terms,
    bool already_sorted_by_magnitude = false);

void sort_frontier_by_magnitude(Frontier& frontier);

// Sorts the frontier by increasing magnitude, with the Pauli key as a
// deterministic tie-breaker, and identifies the largest removable prefix
// whose L1 mass does not exceed cutoff_fraction * total_l1_mass.
L1MassPartition partition_l1_mass(Frontier& frontier, double cutoff_fraction);

std::size_t truncate_l1_mass(Frontier& frontier, double cutoff_fraction);

// Remove every term with |coefficient| < minimum_magnitude, then retain at
// most maximum_support of the remaining largest-magnitude terms. Equal
// magnitudes are ordered by the Pauli key, as in the L1 truncation.
std::size_t truncate_support_budget(
    Frontier& frontier,
    std::size_t maximum_support,
    double minimum_magnitude);

// Compatibility helper for querying the support-budget K without mutation.
std::size_t support_budget_size(
    const Frontier& frontier,
    std::size_t maximum_support,
    double minimum_magnitude);

}  // namespace pauli_bench
