#pragma once

#include "pauli_bench/types.hpp"

namespace pauli_bench {

BfsDiagnostics run_bfs_l1_truncated(const Circuit& circuit, int rz_interval = 4);

}  // namespace pauli_bench
