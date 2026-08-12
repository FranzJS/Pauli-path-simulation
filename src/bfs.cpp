#include "pauli_bench/bfs.hpp"

#include "pauli_bench/kernel.hpp"
#include "pauli_bench/truncation.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace pauli_bench {
namespace {

using Clock = std::chrono::steady_clock;

double terminal_expectation(const Frontier& frontier) {
    double result = 0.0;
    for (const auto& term : frontier) {
        if (term.pauli.x == 0) {
            result += term.coefficient;
        }
    }
    return result;
}

void update_common_peaks(BfsDiagnostics& diagnostics, const Frontier& frontier) {
    diagnostics.peak_support_terms = std::max<std::uint64_t>(
        diagnostics.peak_support_terms,
        frontier.size());
    diagnostics.peak_vector_capacity_terms = std::max<std::uint64_t>(
        diagnostics.peak_vector_capacity_terms,
        frontier.capacity());
}

BfsDiagnostics run_bfs_impl(
    const Circuit& circuit,
    const KStrategyConfig& k_strategy,
    int rz_interval,
    std::vector<std::size_t>* retained_schedule) {
    if (rz_interval <= 0) {
        throw std::invalid_argument("RZ truncation interval must be positive");
    }
    validate_k_strategy(k_strategy);

    BfsDiagnostics diagnostics;
    PauliKernel kernel;
    Frontier frontier;
    frontier.reserve(1);
    frontier.push_back({circuit.observable, 1.0});

    int rz_since_truncation = 0;
    const auto start = Clock::now();
    update_common_peaks(diagnostics, frontier);

    const auto truncate_frontier = [&] {
        diagnostics.peak_pre_truncation_terms = std::max<std::uint64_t>(
            diagnostics.peak_pre_truncation_terms,
            frontier.size());
        const auto decision = determine_k(
            frontier, k_strategy, diagnostics.truncation_events);
        retain_largest_k(
            frontier,
            decision.retained_terms,
            decision.frontier_sorted_by_magnitude);
        if (retained_schedule != nullptr) {
            retained_schedule->push_back(frontier.size());
        }
        diagnostics.peak_post_truncation_terms = std::max<std::uint64_t>(
            diagnostics.peak_post_truncation_terms,
            frontier.size());
        ++diagnostics.truncation_events;
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
            truncate_frontier();
            rz_since_truncation = 0;
        }
    }

    if (rz_since_truncation != 0) {
        truncate_frontier();
    }

    if (k_strategy.strategy == KStrategy::Schedule &&
        diagnostics.truncation_events != k_strategy.schedule.size()) {
        throw std::invalid_argument("K schedule has too many events");
    }

    diagnostics.estimate = terminal_expectation(frontier);
    diagnostics.error = std::isfinite(circuit.reference)
                            ? std::abs(diagnostics.estimate - circuit.reference)
                            : std::numeric_limits<double>::quiet_NaN();
    diagnostics.runtime_seconds =
        std::chrono::duration<double>(Clock::now() - start).count();
    return diagnostics;
}

}  // namespace

BfsDiagnostics run_bfs_truncated(
    const Circuit& circuit,
    const KStrategyConfig& k_strategy,
    int rz_interval,
    std::vector<std::size_t>* retained_schedule) {
    return run_bfs_impl(
        circuit, k_strategy, rz_interval, retained_schedule);
}

BfsDiagnostics run_bfs_l1_truncated(
    const Circuit& circuit,
    int rz_interval,
    std::vector<std::size_t>* retained_schedule) {
    KStrategyConfig config;
    config.strategy = KStrategy::L1Mass;
    config.l1_cutoff = circuit.l1_cutoff;
    return run_bfs_truncated(
        circuit, config, rz_interval, retained_schedule);
}

BfsDiagnostics run_bfs_support_truncated(
    const Circuit& circuit,
    std::size_t maximum_support,
    double minimum_magnitude,
    int rz_interval,
    std::vector<std::size_t>* retained_schedule) {
    KStrategyConfig config;
    config.strategy = KStrategy::SupportBudget;
    config.maximum_support = maximum_support;
    config.minimum_magnitude = minimum_magnitude;
    return run_bfs_truncated(
        circuit, config, rz_interval, retained_schedule);
}

}  // namespace pauli_bench
