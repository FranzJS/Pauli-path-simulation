#pragma once
#include "v1.hpp"
namespace pbv1 {
enum class ResamplingMethod { OrdinaryMultinomial, ResidualMultinomial, ResidualDependentRounding };
struct ResamplingDiagnostics {
 Result result;
 std::uint64_t peak_post_terms{0};
 double average_deterministic_fraction{0};
 std::uint64_t resampling_events{0};
};
const char* resampling_method_name(ResamplingMethod);
ResamplingDiagnostics run_resampling_bfs(const Circuit&, double retained_memory_mb, ResamplingMethod, std::uint64_t seed=20260715);
}
