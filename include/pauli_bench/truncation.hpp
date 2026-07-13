#pragma once
#include "v1.hpp"
namespace pbv1 {
Result run_bfs_l1_truncated(const Circuit&, double budget_s, double discarded_l1_fraction=1e-4);
}
