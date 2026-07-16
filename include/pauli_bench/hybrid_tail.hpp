#pragma once
#include "v1.hpp"
namespace pbv1 {
enum class TailWeighting { OriginalCoefficient, HorvitzThompson };
struct HybridTailDiagnostics {
 Result result;
 std::uint64_t peak_post_terms{0};
 std::uint64_t max_heavy_terms{0};
 std::uint64_t max_sampled_tail_terms{0};
 std::uint64_t truncation_events{0};
 double average_retained_tail_l1_fraction{0};
};
HybridTailDiagnostics run_hybrid_l1_tail(
 const Circuit&, double budget_s, double cutoff, double tail_ratio,
 std::uint64_t seed=20260716, int nonclifford_interval=4,
 TailWeighting weighting=TailWeighting::OriginalCoefficient,
 std::uint64_t support_cap=1500000);
}
