#include "pauli_bench/circuits.hpp"

#include <numbers>
#include <random>
#include <string>

namespace pauli_bench {

Circuit make_clifford_t_brickwork(
    int qubits,
    int layers,
    double t_density,
    std::uint64_t seed,
    double depolarizing_probability) {
    Circuit circuit;
    circuit.family = depolarizing_probability > 0.0 ? "clifford_t_depol" : "clifford_t";
    circuit.name = (depolarizing_probability > 0.0 ? "ct_depol_n" : "ct_n") +
                   std::to_string(qubits) + "_d" + std::to_string(layers);
    circuit.qubits = qubits;
    circuit.observable.x =
        (std::uint64_t{1} << 7) |
        (std::uint64_t{1} << 9) |
        (std::uint64_t{1} << 10);
    circuit.observable.z = std::uint64_t{1} << 8;
    circuit.l1_cutoff = depolarizing_probability > 0.0 ? 0.034 : 0.00625;
    circuit.reference = depolarizing_probability > 0.0
                            ? 0.0663217307060528
                            : 0.329520762689497;
    circuit.reference_method = depolarizing_probability > 0.0
                                   ? "converged_pauli"
                                   : "statevector";

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::uniform_int_distribution<int> clifford_choice(0, 2);

    for (int layer = 0; layer < layers; ++layer) {
        for (int qubit = 0; qubit < qubits; ++qubit) {
            const int choice = clifford_choice(rng);
            if (choice == 0) {
                circuit.gates.push_back(Gate::h(qubit));
            } else if (choice == 1) {
                circuit.gates.push_back(Gate::s(qubit));
            } else {
                circuit.gates.push_back(Gate::h(qubit));
                circuit.gates.push_back(Gate::s(qubit));
            }

            if (uniform(rng) < t_density) {
                circuit.gates.push_back(Gate::rz(qubit, std::numbers::pi / 4.0));
            }
        }

        for (int qubit = layer & 1; qubit + 1 < qubits; qubit += 2) {
            circuit.gates.push_back(Gate::cnot(qubit, qubit + 1));
        }

        if (depolarizing_probability > 0.0) {
            for (int qubit = 0; qubit < qubits; ++qubit) {
                circuit.gates.push_back(
                    Gate::depolarizing(qubit, depolarizing_probability));
            }
        }
    }

    return circuit;
}

Circuit make_nonintegrable_ising(
    int qubits,
    int steps,
    double dt,
    double coupling,
    double transverse_field,
    double longitudinal_field) {
    Circuit circuit;
    circuit.family = "ising";
    circuit.name = "ising_n" + std::to_string(qubits) + "_s" + std::to_string(steps);
    circuit.qubits = qubits;
    circuit.observable.z =
        (std::uint64_t{1} << 9) |
        (std::uint64_t{1} << 10);
    circuit.l1_cutoff = 0.00005;
    circuit.reference = 0.714108793456576;
    circuit.reference_method = "statevector";

    for (int step = 0; step < steps; ++step) {
        for (int qubit = 0; qubit + 1 < qubits; ++qubit) {
            circuit.gates.push_back(Gate::cnot(qubit, qubit + 1));
            circuit.gates.push_back(Gate::rz(qubit + 1, 2.0 * coupling * dt));
            circuit.gates.push_back(Gate::cnot(qubit, qubit + 1));
        }

        for (int qubit = 0; qubit < qubits; ++qubit) {
            circuit.gates.push_back(Gate::h(qubit));
            circuit.gates.push_back(Gate::rz(qubit, 2.0 * transverse_field * dt));
            circuit.gates.push_back(Gate::h(qubit));
            circuit.gates.push_back(Gate::rz(qubit, 2.0 * longitudinal_field * dt));
        }
    }

    return circuit;
}

std::vector<Circuit> benchmark_circuits() {
    return {
        make_clifford_t_brickwork(20, 12, 0.70, 20260715, 0.0),
        make_nonintegrable_ising(20, 12, 0.12, 1.0, 0.91, 0.37),
        make_clifford_t_brickwork(20, 12, 0.70, 20260715, 0.05),
    };
}

}  // namespace pauli_bench
