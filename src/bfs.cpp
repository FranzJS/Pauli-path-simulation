#include "pauli_bench/bfs.hpp"

#include "pauli_bench/kernel.hpp"
#include "pauli_bench/truncation.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>

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

}  // namespace

BfsDiagnostics run_bfs_l1_truncated(
    const Circuit& circuit,
    int rz_interval,
    std::vector<std::size_t>* retained_schedule) {
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
            if (retained_schedule != nullptr) {
                retained_schedule->push_back(frontier.size());
            }
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
        if (retained_schedule != nullptr) {
            retained_schedule->push_back(frontier.size());
        }
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
