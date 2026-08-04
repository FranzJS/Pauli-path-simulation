#include "pauli_bench/bfs.hpp"
#include "pauli_bench/circuits.hpp"

#include <cmath>
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

}  // namespace

int main(int argc, char** argv) {
    if (argc < 6) {
        throw std::invalid_argument(
            "usage: pauli_benchmark QUBITS LAYERS "
            "clifford_t|ising|clifford_t_depol "
            "(l1 L1_CUTOFF | support MAX_SUPPORT MIN_MAGNITUDE | "
            "schedule K_SCHEDULE)");
    }

    const int qubits = std::stoi(argv[1]);
    const int layers = std::stoi(argv[2]);
    const std::string model = argv[3];
    const std::string strategy = argv[4];
    pauli_bench::KStrategyConfig k_strategy;
    if (strategy == "l1") {
        if (argc != 6) {
            throw std::invalid_argument("l1 strategy requires L1_CUTOFF");
        }
        k_strategy.strategy = pauli_bench::KStrategy::L1Mass;
        k_strategy.l1_cutoff = std::stod(argv[5]);
    } else if (strategy == "support") {
        if (argc != 7) {
            throw std::invalid_argument(
                "support strategy requires MAX_SUPPORT and MIN_MAGNITUDE");
        }
        k_strategy.strategy = pauli_bench::KStrategy::SupportBudget;
        k_strategy.maximum_support = parse_maximum_support(argv[5]);
        k_strategy.minimum_magnitude = parse_minimum_magnitude(argv[6]);
    } else if (strategy == "schedule") {
        if (argc != 6) {
            throw std::invalid_argument("schedule strategy requires K_SCHEDULE");
        }
        k_strategy.strategy = pauli_bench::KStrategy::Schedule;
        k_strategy.schedule = parse_schedule(argv[5]);
    } else {
        throw std::invalid_argument("strategy must be l1, support, or schedule");
    }
    pauli_bench::validate_k_strategy(k_strategy);
    const auto circuit = pauli_bench::make_configured_circuit(
        qubits,
        layers,
        model,
        strategy == "l1" ? k_strategy.l1_cutoff : 0.0);
    const auto result = pauli_bench::run_bfs_truncated(
        circuit, k_strategy, 4);
    const std::string method = strategy == "l1"
                                   ? "bfs_l1"
                                   : strategy == "support"
                                         ? "bfs_support"
                                         : "bfs_schedule";

    std::cout
        << circuit.family << ','
        << circuit.name << ','
        << method << ','
        << strategy << ','
        << circuit.qubits << ','
        << std::setprecision(17)
        << (strategy == "l1"
                ? k_strategy.l1_cutoff
                : std::numeric_limits<double>::quiet_NaN()) << ','
        << k_strategy.maximum_support << ','
        << k_strategy.minimum_magnitude << ','
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
        << result.peak_merge_capacity_slots
        << '\n';
}
