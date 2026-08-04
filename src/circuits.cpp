#include "pauli_bench/circuits.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>

namespace pauli_bench {

Circuit make_clifford_t_brickwork(
    int qubits,
    int layers,
    double t_density,
    std::uint64_t seed,
    double depolarizing_probability,
    double l1_cutoff) {
    if (qubits < 4 || qubits > 64) {
        throw std::invalid_argument(
            "Clifford+T circuits require between 4 and 64 qubits");
    }
    if (layers <= 0) {
        throw std::invalid_argument("circuit layers must be positive");
    }
    Circuit circuit;
    circuit.family = depolarizing_probability > 0.0 ? "clifford_t_depol" : "clifford_t";
    circuit.name = (depolarizing_probability > 0.0 ? "ct_depol_n" : "ct_n") +
                   std::to_string(qubits) + "_d" + std::to_string(layers);
    circuit.qubits = qubits;
    const int observable_start = std::clamp(qubits / 2 - 3, 0, qubits - 4);
    circuit.observable.x =
        (std::uint64_t{1} << observable_start) |
        (std::uint64_t{1} << (observable_start + 2)) |
        (std::uint64_t{1} << (observable_start + 3));
    circuit.observable.z = std::uint64_t{1} << (observable_start + 1);
    circuit.l1_cutoff = l1_cutoff >= 0.0
                            ? l1_cutoff
                            : (depolarizing_probability > 0.0 ? 0.034 : 0.00625);
    const bool has_fixed_reference =
        qubits == 20 && layers == 12 && t_density == 0.70 &&
        seed == 20260715 &&
        (depolarizing_probability == 0.0 || depolarizing_probability == 0.05);
    if (has_fixed_reference) {
        circuit.reference = depolarizing_probability > 0.0
                                ? 0.0663217307060528
                                : 0.329520762689497;
        circuit.reference_method = depolarizing_probability > 0.0
                                       ? "converged_pauli"
                                       : "statevector";
    } else {
        circuit.reference = std::numeric_limits<double>::quiet_NaN();
        circuit.reference_method = "unavailable";
    }

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
    double longitudinal_field,
    double l1_cutoff) {
    if (qubits < 2 || qubits > 64) {
        throw std::invalid_argument(
            "Ising circuits require between 2 and 64 qubits");
    }
    if (steps <= 0) {
        throw std::invalid_argument("Trotter steps must be positive");
    }
    Circuit circuit;
    circuit.family = "ising";
    circuit.name = "ising_n" + std::to_string(qubits) + "_s" + std::to_string(steps);
    circuit.qubits = qubits;
    const int observable_start = std::clamp(qubits / 2 - 1, 0, qubits - 2);
    circuit.observable.z =
        (std::uint64_t{1} << observable_start) |
        (std::uint64_t{1} << (observable_start + 1));
    circuit.l1_cutoff = l1_cutoff >= 0.0 ? l1_cutoff : 0.00005;
    const bool has_fixed_reference =
        qubits == 20 && steps == 12 && dt == 0.12 && coupling == 1.0 &&
        transverse_field == 0.91 && longitudinal_field == 0.37;
    if (has_fixed_reference) {
        circuit.reference = 0.714108793456576;
        circuit.reference_method = "statevector";
    } else {
        circuit.reference = std::numeric_limits<double>::quiet_NaN();
        circuit.reference_method = "unavailable";
    }

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

Circuit make_configured_circuit(
    int qubits,
    int layers,
    const std::string& model,
    double l1_cutoff) {
    if (!std::isfinite(l1_cutoff) || l1_cutoff < 0.0 || l1_cutoff > 1.0) {
        throw std::invalid_argument("L1 cutoff must be finite and in [0, 1]");
    }
    if (model == "clifford_t") {
        return make_clifford_t_brickwork(
            qubits, layers, 0.70, 20260715, 0.0, l1_cutoff);
    }
    if (model == "clifford_t_depol") {
        return make_clifford_t_brickwork(
            qubits, layers, 0.70, 20260715, 0.05, l1_cutoff);
    }
    if (model == "ising") {
        return make_nonintegrable_ising(
            qubits, layers, 0.12, 1.0, 0.91, 0.37, l1_cutoff);
    }
    throw std::invalid_argument(
        "model must be clifford_t, ising, or clifford_t_depol");
}

std::vector<Circuit> benchmark_circuits() {
    return {
        make_configured_circuit(20, 12, "clifford_t", 0.00625),
        make_configured_circuit(20, 12, "ising", 0.00005),
        make_configured_circuit(20, 12, "clifford_t_depol", 0.034),
    };
}

}  // namespace pauli_bench
