#include "pauli_bench/bfs.hpp"
#include "pauli_bench/circuits.hpp"
#include "pauli_bench/hybrid.hpp"

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::size_t parse_maximum_support(const std::string& encoded) {
    if (encoded.empty() || encoded.front() == '-') {
        throw std::invalid_argument("MAX_SUPPORT must be a nonnegative integer");
    }
    return static_cast<std::size_t>(std::stoull(encoded));
}

double parse_minimum_magnitude(const std::string& encoded) {
    const double value = std::stod(encoded);
    if (!std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument(
            "MIN_MAGNITUDE must be finite and nonnegative");
    }
    return value;
}

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
    if (schedule.empty()) {
        throw std::invalid_argument("K_SCHEDULE must not be empty");
    }
    return schedule;
}

void print_result(
    const pauli_bench::Circuit& circuit,
    const std::string& method,
    const std::string& truncation_strategy,
    double l1_cutoff,
    std::uint64_t maximum_support,
    double minimum_magnitude,
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
        << truncation_strategy << ','
        << std::setprecision(17) << tail_ratio << ','
        << seed << ','
        << circuit.qubits << ','
        << l1_cutoff << ','
        << maximum_support << ','
        << minimum_magnitude << ','
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
    if (argc < 7) {
        throw std::invalid_argument(
            "usage: pauli_magnitude_pps_benchmark QUBITS LAYERS "
            "clifford_t|ising|clifford_t_depol "
            "(l1 L1_CUTOFF | support MAX_SUPPORT MIN_MAGNITUDE | "
            "schedule K_SCHEDULE) bfs|schedule|pps_ht [SEED]");
    }

    const int qubits = std::stoi(argv[1]);
    const int layers = std::stoi(argv[2]);
    const std::string model = argv[3];
    const std::string strategy = argv[4];
    pauli_bench::KStrategyConfig k_strategy;
    int method_index{};
    if (strategy == "l1") {
        k_strategy.strategy = pauli_bench::KStrategy::L1Mass;
        k_strategy.l1_cutoff = std::stod(argv[5]);
        method_index = 6;
    } else if (strategy == "support") {
        if (argc < 8) {
            throw std::invalid_argument(
                "support strategy requires MAX_SUPPORT and MIN_MAGNITUDE");
        }
        k_strategy.strategy = pauli_bench::KStrategy::SupportBudget;
        k_strategy.maximum_support = parse_maximum_support(argv[5]);
        k_strategy.minimum_magnitude = parse_minimum_magnitude(argv[6]);
        method_index = 7;
    } else if (strategy == "schedule") {
        k_strategy.strategy = pauli_bench::KStrategy::Schedule;
        k_strategy.schedule = parse_schedule(argv[5]);
        method_index = 6;
    } else {
        throw std::invalid_argument("strategy must be l1, support, or schedule");
    }
    pauli_bench::validate_k_strategy(k_strategy);
    const std::string method = argv[method_index];
    const auto circuit = pauli_bench::make_configured_circuit(
        qubits,
        layers,
        model,
        strategy == "l1" ? k_strategy.l1_cutoff : 0.0);

    if (method == "bfs") {
        const auto result = pauli_bench::run_bfs_truncated(
            circuit, k_strategy, 4);
        const std::string bfs_method = strategy == "l1"
                                           ? "bfs_l1"
                                           : strategy == "support"
                                                 ? "bfs_support"
                                                 : "bfs_schedule";
        print_result(
            circuit,
            bfs_method,
            strategy,
            strategy == "l1"
                ? k_strategy.l1_cutoff
                : std::numeric_limits<double>::quiet_NaN(),
            k_strategy.maximum_support,
            k_strategy.minimum_magnitude,
            0.0, 0, result, 0, 0, 0, 1.0, 0.0);
        return 0;
    }
    if (method == "schedule") {
        std::vector<std::size_t> retained_schedule;
        pauli_bench::run_bfs_truncated(
            circuit, k_strategy, 4, &retained_schedule);
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
        if (argc <= method_index + 1) {
            throw std::invalid_argument("pps_ht requires SEED");
        }
        const auto seed = static_cast<std::uint64_t>(
            std::stoull(argv[method_index + 1]));
        if (argc != method_index + 2) {
            throw std::invalid_argument(
                "pps_ht accepts exactly one SEED after the method");
        }
        const auto result = pauli_bench::run_optimal_pps_ht(
            circuit, seed, k_strategy, 4);
        print_result(
            circuit,
            "l1_optimal_pps_ht",
            strategy,
            strategy == "l1"
                ? k_strategy.l1_cutoff
                : std::numeric_limits<double>::quiet_NaN(),
            k_strategy.maximum_support,
            k_strategy.minimum_magnitude,
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
