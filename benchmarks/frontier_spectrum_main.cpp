#include "pauli_bench/circuits.hpp"
#include "pauli_bench/kernel.hpp"
#include "pauli_bench/reference.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

pauli_bench::Frontier first_exact_truncation_frontier_at_least(
    const pauli_bench::Circuit& circuit,
    int rz_interval,
    std::size_t minimum_terms,
    std::size_t& event_index) {
    pauli_bench::PauliKernel kernel;
    pauli_bench::Frontier frontier{{circuit.observable, 1.0}};
    int rz_since_event = 0;
    std::size_t event = 0;
    event_index = 0;

    for (auto gate = circuit.gates.rbegin(); gate != circuit.gates.rend(); ++gate) {
        if (gate->kind == pauli_bench::GateKind::RZ) {
            kernel.apply_rz_and_merge(frontier, *gate);
            ++rz_since_event;
        } else {
            kernel.apply_deterministic_in_place(frontier, *gate);
        }
        if (rz_since_event == rz_interval) {
            if (frontier.size() >= minimum_terms) {
                event_index = event;
                return frontier;
            }
            ++event;
            rz_since_event = 0;
        }
    }
    if (rz_since_event != 0) {
        event_index = event;
    }
    return frontier;
}

bool same_plateau(double lhs, double rhs) {
    const double scale = std::max({std::abs(lhs), std::abs(rhs), 1e-300});
    return std::abs(lhs - rhs) <= 1e-11 * scale;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        throw std::invalid_argument(
            "usage: pauli_frontier_spectrum QUBITS LAYERS clifford_t MIN_TERMS");
    }
    const int qubits = std::stoi(argv[1]);
    const int layers = std::stoi(argv[2]);
    const std::string model = argv[3];
    const auto minimum_terms = static_cast<std::size_t>(std::stoull(argv[4]));
    if (model != "clifford_t") {
        throw std::invalid_argument(
            "frontier spectrum example currently requires clifford_t");
    }

    auto circuit = pauli_bench::make_configured_circuit(
        qubits, layers, model, 0.0);
    const double reference = pauli_bench::statevector_expectation(circuit);
    std::size_t event_index = 0;
    auto frontier = first_exact_truncation_frontier_at_least(
        circuit, 4, minimum_terms, event_index);

    std::vector<double> magnitudes;
    magnitudes.reserve(frontier.size());
    for (const auto& term : frontier) {
        magnitudes.push_back(std::abs(term.coefficient));
    }
    std::sort(magnitudes.begin(), magnitudes.end(), std::greater<double>());

    std::cout << std::setprecision(17);
    std::cout << "reference," << reference << '\n';
    std::cout << "event," << event_index << '\n';
    std::cout << "frontier_terms," << magnitudes.size() << '\n';
    std::cout << "rank_start,rank_end,magnitude,count\n";
    std::size_t start = 0;
    while (start < magnitudes.size()) {
        std::size_t end = start + 1;
        while (end < magnitudes.size() &&
               same_plateau(magnitudes[start], magnitudes[end])) {
            ++end;
        }
        std::cout << start + 1 << ',' << end << ',' << magnitudes[start]
                  << ',' << end - start << '\n';
        start = end;
    }
}
