#pragma once
#include "v1.hpp"
namespace pbv1 {
Result run_terminal_aware_monte_carlo(const Circuit&, double budget_s, std::uint64_t seed=20260715);
}
