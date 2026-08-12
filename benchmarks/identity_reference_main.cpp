#include "pauli_bench/circuits.hpp"

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    if (argc != 5) {
        throw std::invalid_argument(
            "usage: pauli_identity_reference QUBITS W_DEPTH "
            "clifford_t|ising CIRCUIT_SEED");
    }
    const int qubits = std::stoi(argv[1]);
    const int prefix_layers = std::stoi(argv[2]);
    const std::string model = argv[3];
    const auto circuit_seed = static_cast<std::uint64_t>(
        std::stoull(argv[4]));
    const auto reference = pauli_bench::compute_prefixed_identity_reference(
        qubits, prefix_layers, model, circuit_seed);
    std::cout << std::setprecision(17) << reference.value << ','
              << reference.lightcone_qubits << '\n';
}
