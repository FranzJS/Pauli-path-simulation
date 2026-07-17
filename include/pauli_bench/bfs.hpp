#pragma once

#include "pauli_bench/types.hpp"

#include <cstddef>
#include <vector>

namespace pauli_bench {

BfsDiagnostics run_bfs_l1_truncated(
    const Circuit& circuit,
    int rz_interval = 4,
    std::vector<std::size_t>* retained_schedule = nullptr);

}  // namespace pauli_bench
