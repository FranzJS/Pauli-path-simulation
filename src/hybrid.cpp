#include "pauli_bench/hybrid.hpp"

#include "pauli_bench/kernel.hpp"
#include "pauli_bench/truncation.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace pauli_bench {
namespace {

using Clock = std::chrono::steady_clock;

struct PpsSample {
    std::vector<bool> selected;
    std::vector<long double> inclusion_probability;
    std::size_t saturated_terms{};
};

template <typename WeightAt>
PpsSample sample_fixed_size_pps(
    const std::size_t population_terms,
    std::size_t sample_terms,
    std::mt19937_64& rng,
    WeightAt weight_at,
    bool weights_sorted_ascending = false) {
    PpsSample sample{
        std::vector<bool>(population_terms, false),
        std::vector<long double>(population_terms, 0.0L)};
    if (population_terms == 0 || sample_terms == 0) {
        return sample;
    }
    sample_terms = std::min(sample_terms, population_terms);
    if (sample_terms == population_terms) {
        std::fill(sample.selected.begin(), sample.selected.end(), true);
        std::fill(
            sample.inclusion_probability.begin(),
            sample.inclusion_probability.end(),
            1.0L);
        sample.saturated_terms = population_terms;
        return sample;
    }

    long double remaining_mass = 0.0L;
    for (std::size_t index = 0; index < population_terms; ++index) {
        const long double weight = weight_at(index);
        if (!std::isfinite(weight) || weight < 0.0L) {
            throw std::invalid_argument(
                "PPS sampling weights must be finite and nonnegative");
        }
        remaining_mass += weight;
    }

    std::vector<std::size_t> sorted_by_weight;
    if (!weights_sorted_ascending) {
        sorted_by_weight.resize(population_terms);
        std::iota(sorted_by_weight.begin(), sorted_by_weight.end(), 0);
        std::sort(
            sorted_by_weight.begin(),
            sorted_by_weight.end(),
            [&](const std::size_t lhs, const std::size_t rhs) {
                const long double lhs_weight = weight_at(lhs);
                const long double rhs_weight = weight_at(rhs);
                if (lhs_weight != rhs_weight) {
                    return lhs_weight < rhs_weight;
                }
                return lhs < rhs;
            });
    }
    const auto ranked_index = [&](const std::size_t position) {
        return weights_sorted_ascending
                   ? position
                   : sorted_by_weight[position];
    };

    std::vector<std::size_t> order(population_terms);
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), rng);

    if (remaining_mass <= 0.0L) {
        const long double probability =
            static_cast<long double>(sample_terms) /
            static_cast<long double>(population_terms);
        std::fill(
            sample.inclusion_probability.begin(),
            sample.inclusion_probability.end(),
            probability);
        for (std::size_t index = 0; index < sample_terms; ++index) {
            sample.selected[order[index]] = true;
        }
        return sample;
    }

    // Find tau such that pi_i = min(1, |c_i| / tau) and sum_i pi_i = M.
    // The population is already sorted by increasing |c_i|.
    std::size_t saturated = 0;
    long double tau = 0.0L;
    while (true) {
        const std::size_t unsaturated = population_terms - saturated;
        const std::size_t remaining_slots = sample_terms - saturated;
        if (unsaturated == 0 || remaining_slots == 0) {
            break;
        }
        tau = remaining_mass / static_cast<long double>(remaining_slots);
        const auto largest_unsaturated = unsaturated - 1;
        const long double largest_weight =
            weight_at(ranked_index(largest_unsaturated));
        if (largest_weight <= tau * (1.0L + 1e-18L)) {
            break;
        }
        remaining_mass -= largest_weight;
        ++saturated;
    }
    sample.saturated_terms = saturated;

    long double probability_sum = 0.0L;
    const std::size_t unsaturated = population_terms - saturated;
    const std::size_t remaining_slots = sample_terms - saturated;
    for (std::size_t index = unsaturated;
         index < population_terms;
         ++index) {
        sample.inclusion_probability[ranked_index(index)] = 1.0L;
        probability_sum += 1.0L;
    }
    if (remaining_slots != 0) {
        if (remaining_mass > 0.0L) {
            tau = remaining_mass / static_cast<long double>(remaining_slots);
            for (std::size_t position = 0;
                 position < unsaturated;
                 ++position) {
                const auto index = ranked_index(position);
                sample.inclusion_probability[index] = std::min<long double>(
                    1.0L,
                    weight_at(index) / tau);
                probability_sum += sample.inclusion_probability[index];
            }
        } else {
            const long double probability =
                static_cast<long double>(remaining_slots) /
                static_cast<long double>(unsaturated);
            for (std::size_t position = 0;
                 position < unsaturated;
                 ++position) {
                const auto index = ranked_index(position);
                sample.inclusion_probability[index] = probability;
                probability_sum += probability;
            }
        }
    }

    // Correct only floating-point drift so systematic sampling sees exactly M
    // units of probability mass.
    long double drift =
        static_cast<long double>(sample_terms) - probability_sum;
    for (const auto index : order) {
        if (std::abs(drift) <= 1e-14L) {
            break;
        }
        const long double room = drift > 0.0L
                                     ? 1.0L - sample.inclusion_probability[index]
                                     : sample.inclusion_probability[index];
        const long double adjustment = std::copysign(
            std::min(std::abs(drift), room),
            drift);
        sample.inclusion_probability[index] += adjustment;
        drift -= adjustment;
    }
    probability_sum = std::accumulate(
        sample.inclusion_probability.begin(),
        sample.inclusion_probability.end(),
        0.0L);
    drift = static_cast<long double>(sample_terms) - probability_sum;
    if (drift != 0.0L) {
        for (const auto index : order) {
            const long double adjusted =
                sample.inclusion_probability[index] + drift;
            if (adjusted >= 0.0L && adjusted <= 1.0L) {
                sample.inclusion_probability[index] = adjusted;
                drift = 0.0L;
                break;
            }
        }
    }
    if (drift != 0.0L) {
        throw std::runtime_error("failed to normalize PPS inclusion probabilities");
    }

    long double point = std::generate_canonical<
        long double,
        std::numeric_limits<long double>::digits>(rng);
    long double cumulative = 0.0L;
    std::size_t selected_count = 0;
    for (std::size_t position = 0; position < order.size(); ++position) {
        const auto index = order[position];
        const long double next = position + 1 == order.size()
                                     ? static_cast<long double>(sample_terms)
                                     : cumulative +
                                           sample.inclusion_probability[index];
        if (selected_count < sample_terms && point < next) {
            sample.selected[index] = true;
            ++selected_count;
            point += 1.0L;
        }
        cumulative = next;
    }

    if (selected_count != sample_terms) {
        throw std::runtime_error("systematic PPS produced the wrong sample size");
    }
    return sample;
}

double terminal_expectation(const Frontier& frontier) {
    double result = 0.0;
    for (const auto& term : frontier) {
        if (term.pauli.x == 0) {
            result += term.coefficient;
        }
    }
    return result;
}

void update_common_peaks(
    HybridDiagnostics& diagnostics,
    const Frontier& frontier) {
    diagnostics.peak_support_terms = std::max<std::uint64_t>(
        diagnostics.peak_support_terms,
        frontier.size());
    diagnostics.peak_vector_capacity_terms = std::max<std::uint64_t>(
        diagnostics.peak_vector_capacity_terms,
        frontier.capacity());
}

}  // namespace

template <typename WeightAt>
HybridTruncationStats truncate_weighted_pps_ht_impl(
    Frontier& frontier,
    std::size_t target_total_terms,
    std::mt19937_64& rng,
    WeightAt weight_at,
    bool weights_sorted_ascending = false) {
    HybridTruncationStats stats;
    const std::size_t sample_terms =
        std::min(target_total_terms, frontier.size());
    const auto sample = sample_fixed_size_pps(
        frontier.size(),
        sample_terms,
        rng,
        weight_at,
        weights_sorted_ascending);

    stats.heavy_terms = sample.saturated_terms;
    stats.tail_terms = frontier.size() - stats.heavy_terms;
    stats.sampled_tail_terms = sample_terms - stats.heavy_terms;
    for (std::size_t index = 0; index < frontier.size(); ++index) {
        if (sample.inclusion_probability[index] < 1.0L) {
            stats.tail_l1_mass += std::abs(frontier[index].coefficient);
        }
    }

    std::size_t write = 0;
    for (std::size_t index = 0; index < frontier.size(); ++index) {
        if (!sample.selected[index]) {
            continue;
        }
        const long double probability = sample.inclusion_probability[index];
        if (probability <= 0.0L) {
            throw std::runtime_error(
                "selected PPS coordinate has zero inclusion probability");
        }
        const double multiplier = 1.0 / static_cast<double>(probability);
        stats.max_importance_multiplier = std::max(
            stats.max_importance_multiplier,
            multiplier);
        frontier[index].coefficient *= multiplier;
        frontier[write++] = std::move(frontier[index]);
    }
    frontier.resize(write);
    if (frontier.size() != sample_terms) {
        throw std::runtime_error("optimal PPS truncation violated budget K");
    }
    return stats;
}

HybridTruncationStats truncate_l1_optimal_pps_ht_impl(
    Frontier& frontier,
    std::size_t target_total_terms,
    std::mt19937_64& rng) {
    sort_frontier_by_magnitude(frontier);
    return truncate_weighted_pps_ht_impl(
        frontier,
        target_total_terms,
        rng,
        [&](const std::size_t index) {
            return static_cast<long double>(
                std::abs(frontier[index].coefficient));
        },
        true);
}

HybridDiagnostics run_l1_optimal_pps_ht_impl(
    const Circuit& circuit,
    std::uint64_t seed,
    const std::vector<std::size_t>& retained_schedule,
    int rz_interval) {
    if (rz_interval <= 0) {
        throw std::invalid_argument("RZ truncation interval must be positive");
    }

    HybridDiagnostics diagnostics;
    PauliKernel kernel;
    Frontier frontier;
    frontier.reserve(1);
    frontier.push_back({circuit.observable, 1.0});
    std::mt19937_64 rng(seed);

    int rz_since_truncation = 0;
    std::size_t truncation_index = 0;
    const auto start = Clock::now();
    update_common_peaks(diagnostics, frontier);

    const auto truncate = [&] {
        diagnostics.peak_pre_truncation_terms = std::max<std::uint64_t>(
            diagnostics.peak_pre_truncation_terms,
            frontier.size());
        if (truncation_index >= retained_schedule.size()) {
            throw std::invalid_argument("retained schedule has too few events");
        }
        const auto stats = truncate_l1_optimal_pps_ht_impl(
            frontier,
            retained_schedule[truncation_index],
            rng);
        ++truncation_index;
        diagnostics.max_heavy_terms = std::max<std::uint64_t>(
            diagnostics.max_heavy_terms,
            stats.heavy_terms);
        diagnostics.max_sampled_tail_terms = std::max<std::uint64_t>(
            diagnostics.max_sampled_tail_terms,
            stats.sampled_tail_terms);
        diagnostics.total_sampled_tail_terms += stats.sampled_tail_terms;
        diagnostics.max_importance_multiplier = std::max(
            diagnostics.max_importance_multiplier,
            stats.max_importance_multiplier);
        for (const auto& term : frontier) {
            diagnostics.max_post_truncation_abs_coefficient = std::max(
                diagnostics.max_post_truncation_abs_coefficient,
                std::abs(term.coefficient));
        }
        diagnostics.peak_post_truncation_terms = std::max<std::uint64_t>(
            diagnostics.peak_post_truncation_terms,
            frontier.size());
        ++diagnostics.truncation_events;
        update_common_peaks(diagnostics, frontier);
    };

    for (auto gate = circuit.gates.rbegin(); gate != circuit.gates.rend(); ++gate) {
        if (gate->kind == GateKind::RZ) {
            kernel.apply_rz_and_merge(frontier, *gate);
            diagnostics.peak_merge_entries = std::max<std::uint64_t>(
                diagnostics.peak_merge_entries,
                kernel.merge_entries());
            diagnostics.peak_merge_capacity_slots = std::max<std::uint64_t>(
                diagnostics.peak_merge_capacity_slots,
                kernel.merge_capacity());
            ++rz_since_truncation;
        } else {
            kernel.apply_deterministic_in_place(frontier, *gate);
        }

        update_common_peaks(diagnostics, frontier);
        if (rz_since_truncation == rz_interval) {
            truncate();
            rz_since_truncation = 0;
        }
    }

    if (rz_since_truncation != 0) {
        truncate();
    }

    if (truncation_index != retained_schedule.size()) {
        throw std::invalid_argument("retained schedule has too many events");
    }

    diagnostics.estimate = terminal_expectation(frontier);
    diagnostics.error = std::abs(diagnostics.estimate - circuit.reference);
    diagnostics.runtime_seconds =
        std::chrono::duration<double>(Clock::now() - start).count();
    return diagnostics;
}

HybridTruncationStats truncate_l1_optimal_pps_ht(
    Frontier& frontier,
    std::size_t target_total_terms,
    std::mt19937_64& rng) {
    return truncate_l1_optimal_pps_ht_impl(
        frontier,
        target_total_terms,
        rng);
}

HybridDiagnostics run_l1_optimal_pps_ht(
    const Circuit& circuit,
    std::uint64_t seed,
    const std::vector<std::size_t>& retained_schedule,
    int rz_interval) {
    return run_l1_optimal_pps_ht_impl(
        circuit,
        seed,
        retained_schedule,
        rz_interval);
}

}  // namespace pauli_bench
