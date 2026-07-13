#pragma once
#include "types.hpp"
namespace pb {
Result run_statevector(const Circuit&, double budget_s);
Result run_exact_sparse(const Circuit&, double budget_s, std::uint64_t term_cap=500000);
Result run_dfs(const Circuit&, double budget_s, std::uint64_t path_cap=200000000);
Result run_monte_carlo(const Circuit&, double budget_s, std::uint64_t max_samples=1000000, std::uint64_t seed=7);
}
