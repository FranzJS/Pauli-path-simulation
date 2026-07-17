#include "pauli_bench/bfs.hpp"
#include "pauli_bench/circuits.hpp"
#include "pauli_bench/hybrid.hpp"

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<std::size_t> parse_schedule(const std::string& encoded) {
    std::vector<std::size_t> schedule;
    std::size_t begin = 0;
    while (begin < encoded.size()) {
        const auto end = encoded.find(':', begin);
        schedule.push_back(static_cast<std::size_t>(std::stoull(
            encoded.substr(begin, end - begin))));
        if (end == std::string::npos) {
            break;
        }
        begin = end + 1;
    }
    return schedule;
}

void print_result(
    const pauli_bench::Circuit& circuit,
    const std::string& method,
    double tail_ratio,
    std::uint64_t seed,
    const pauli_bench::BfsDiagnostics& result,
    std::uint64_t max_heavy_terms,
    std::uint64_t max_sampled_tail_terms,
    std::uint64_t total_sampled_tail_terms,
    double max_importance_multiplier,
    double max_post_truncation_abs_coefficient) {
    std::cout
        << circuit.family << ','
        << circuit.name << ','
        << method << ','
        << std::setprecision(17) << tail_ratio << ','
        << seed << ','
        << circuit.qubits << ','
        << circuit.l1_cutoff << ','
        << result.estimate << ','
        << circuit.reference << ','
        << circuit.reference_method << ','
        << result.error << ','
        << result.runtime_seconds << ','
        << result.peak_support_terms << ','
        << result.peak_pre_truncation_terms << ','
        << result.peak_post_truncation_terms << ','
        << result.truncation_events << ','
        << result.peak_vector_capacity_terms << ','
        << result.peak_merge_entries << ','
        << result.peak_merge_capacity_slots << ','
        << max_heavy_terms << ','
        << max_sampled_tail_terms << ','
        << total_sampled_tail_terms << ','
        << max_importance_multiplier << ','
        << max_post_truncation_abs_coefficient
        << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        throw std::invalid_argument(
            "usage: pauli_magnitude_pps_benchmark CASE "
            "bfs|schedule|pps_ht [SEED] [K_SCHEDULE]");
    }

    const int case_index = std::stoi(argv[1]);
    const std::string method = argv[2];
    const auto circuits = pauli_bench::benchmark_circuits();
    if (case_index < 0 || case_index >= static_cast<int>(circuits.size())) {
        throw std::out_of_range("case index must be 0, 1, or 2");
    }
    const auto& circuit = circuits[case_index];

    if (method == "bfs") {
        const auto result = pauli_bench::run_bfs_l1_truncated(circuit, 4);
        print_result(circuit, "bfs_l1", 0.0, 0, result, 0, 0, 0, 1.0, 0.0);
        return 0;
    }
    if (method == "schedule") {
        std::vector<std::size_t> retained_schedule;
        pauli_bench::run_bfs_l1_truncated(circuit, 4, &retained_schedule);
        for (std::size_t index = 0; index < retained_schedule.size(); ++index) {
            if (index != 0) {
                std::cout << ':';
            }
            std::cout << retained_schedule[index];
        }
        std::cout << '\n';
        return 0;
    }
    if (method == "pps_ht") {
        if (argc < 5) {
            throw std::invalid_argument(
                "pps_ht requires SEED and K_SCHEDULE");
        }
        const auto seed = static_cast<std::uint64_t>(std::stoull(argv[3]));
        const auto retained_schedule = parse_schedule(argv[4]);
        const auto result = pauli_bench::run_l1_optimal_pps_ht(
            circuit,
            seed,
            retained_schedule,
            4);
        print_result(
            circuit,
            "l1_optimal_pps_ht",
            1.0,
            seed,
            result,
            result.max_heavy_terms,
            result.max_sampled_tail_terms,
            result.total_sampled_tail_terms,
            result.max_importance_multiplier,
            result.max_post_truncation_abs_coefficient);
        return 0;
    }
    throw std::invalid_argument("method must be bfs, schedule, or pps_ht");
}
