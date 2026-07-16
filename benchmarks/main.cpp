#include "pauli_bench/bfs.hpp"
#include "pauli_bench/circuits.hpp"

#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    const int case_index = argc > 1 ? std::stoi(argv[1]) : 0;
    const auto circuits = pauli_bench::benchmark_circuits();
    if (case_index < 0 || case_index >= static_cast<int>(circuits.size())) {
        throw std::out_of_range("case index must be 0, 1, or 2");
    }

    const auto& circuit = circuits[case_index];
    const auto result = pauli_bench::run_bfs_l1_truncated(circuit, 4);

    std::cout
        << circuit.family << ','
        << circuit.name << ','
        << "bfs_l1" << ','
        << circuit.qubits << ','
        << std::setprecision(17) << circuit.l1_cutoff << ','
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
