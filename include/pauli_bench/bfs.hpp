#pragma once

#include "pauli_bench/truncation.hpp"
#include "pauli_bench/types.hpp"

#include <cstddef>
#include <vector>

namespace pauli_bench {

BfsDiagnostics run_bfs_truncated(
    const Circuit& circuit,
    const KStrategyConfig& k_strategy,
    int rz_interval = 4,
    std::vector<std::size_t>* retained_schedule = nullptr);

BfsDiagnostics run_bfs_l1_truncated(
    const Circuit& circuit,
    int rz_interval = 4,
    std::vector<std::size_t>* retained_schedule = nullptr);

BfsDiagnostics run_bfs_support_truncated(
    const Circuit& circuit,
    std::size_t maximum_support,
    double minimum_magnitude,
    int rz_interval = 4,
    std::vector<std::size_t>* retained_schedule = nullptr);

}  // namespace pauli_bench
