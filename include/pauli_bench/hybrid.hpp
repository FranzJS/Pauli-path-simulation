#pragma once

#include "pauli_bench/truncation.hpp"
#include "pauli_bench/types.hpp"

#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

namespace pauli_bench {

struct HybridTruncationStats {
    std::size_t heavy_terms{};
    std::size_t tail_terms{};
    std::size_t sampled_tail_terms{};
    double tail_l1_mass{};
    double max_importance_multiplier{1.0};
};

struct HybridDiagnostics : BfsDiagnostics {
    std::uint64_t max_heavy_terms{};
    std::uint64_t max_sampled_tail_terms{};
    std::uint64_t total_sampled_tail_terms{};
    double max_importance_multiplier{1.0};
    double max_post_truncation_abs_coefficient{};
};

// Selects exactly K coordinates from the entire frontier. The first-order
// inclusion probabilities minimize conditional squared L2 error:
//     pi_i = min(1, |c_i| / tau),  sum_i pi_i = K.
// Coordinates with pi_i == 1 are retained deterministically; every other
// selected coefficient receives the Horvitz--Thompson correction 1 / pi_i.
HybridTruncationStats truncate_l1_optimal_pps_ht(
    Frontier& frontier,
    std::size_t target_total_terms,
    std::mt19937_64& rng);

HybridDiagnostics run_optimal_pps_ht(
    const Circuit& circuit,
    std::uint64_t seed,
    const KStrategyConfig& k_strategy,
    int rz_interval = 4);

HybridDiagnostics run_l1_optimal_pps_ht(
    const Circuit& circuit,
    std::uint64_t seed,
    const std::vector<std::size_t>& retained_schedule,
    int rz_interval = 4);

// At each truncation event choose
//   K = min(maximum_support, count(|coefficient| >= minimum_magnitude))
// from the randomized pass's current frontier, then apply the same whole-
// frontier PPS/HT sampler used by the scheduled method.
HybridDiagnostics run_l1_optimal_pps_ht_support_budget(
    const Circuit& circuit,
    std::uint64_t seed,
    std::size_t maximum_support,
    double minimum_magnitude,
    int rz_interval = 4);

}  // namespace pauli_bench
