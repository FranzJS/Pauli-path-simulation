#include "pauli_bench/bfs.hpp"

#include "pauli_bench/kernel.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>

namespace pauli_bench {
namespace {

using Clock = std::chrono::steady_clock;

std::size_t truncate_l1_mass(Frontier& frontier, double cutoff_fraction) {
    if (frontier.empty() || cutoff_fraction <= 0.0) {
        return 0;
    }

    double l1_mass = 0.0;
    for (const auto& term : frontier) {
        l1_mass += std::abs(term.coefficient);
    }

    const double removal_budget = cutoff_fraction * l1_mass;
    if (removal_budget <= 0.0) {
        return 0;
    }

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

    double removed_mass = 0.0;
    std::size_t removed_terms = 0;
    while (removed_terms < frontier.size()) {
        const double next_mass = std::abs(frontier[removed_terms].coefficient);
        if (removed_mass + next_mass > removal_budget) {
            break;
        }
        removed_mass += next_mass;
        ++removed_terms;
    }

    if (removed_terms == 0) {
        return 0;
    }

    std::move(
        frontier.begin() + static_cast<std::ptrdiff_t>(removed_terms),
        frontier.end(),
        frontier.begin());
    frontier.resize(frontier.size() - removed_terms);
    return removed_terms;
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

void update_common_peaks(BfsDiagnostics& diagnostics, const Frontier& frontier) {
    diagnostics.peak_support_terms = std::max<std::uint64_t>(
        diagnostics.peak_support_terms,
        frontier.size());
    diagnostics.peak_vector_capacity_terms = std::max<std::uint64_t>(
        diagnostics.peak_vector_capacity_terms,
        frontier.capacity());
}

}  // namespace

BfsDiagnostics run_bfs_l1_truncated(const Circuit& circuit, int rz_interval) {
    BfsDiagnostics diagnostics;
    PauliKernel kernel;
    Frontier frontier;
    frontier.reserve(1);
    frontier.push_back({circuit.observable, 1.0});

    int rz_since_truncation = 0;
    const auto start = Clock::now();

    update_common_peaks(diagnostics, frontier);

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
            diagnostics.peak_pre_truncation_terms = std::max<std::uint64_t>(
                diagnostics.peak_pre_truncation_terms,
                frontier.size());
            truncate_l1_mass(frontier, circuit.l1_cutoff);
            diagnostics.peak_post_truncation_terms = std::max<std::uint64_t>(
                diagnostics.peak_post_truncation_terms,
                frontier.size());
            ++diagnostics.truncation_events;
            rz_since_truncation = 0;
        }
    }

    if (rz_since_truncation != 0) {
        diagnostics.peak_pre_truncation_terms = std::max<std::uint64_t>(
            diagnostics.peak_pre_truncation_terms,
            frontier.size());
        truncate_l1_mass(frontier, circuit.l1_cutoff);
        diagnostics.peak_post_truncation_terms = std::max<std::uint64_t>(
            diagnostics.peak_post_truncation_terms,
            frontier.size());
        ++diagnostics.truncation_events;
    }

    diagnostics.estimate = terminal_expectation(frontier);
    diagnostics.error = std::abs(diagnostics.estimate - circuit.reference);
    diagnostics.runtime_seconds =
        std::chrono::duration<double>(Clock::now() - start).count();
    return diagnostics;
}

}  // namespace pauli_bench
