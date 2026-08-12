#include "pauli_bench/circuits.hpp"
#include "pauli_bench/reference.hpp"
#include "pauli_bench/reference_registry.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>

namespace pauli_bench {
namespace {

std::uint64_t prefix_circuit_seed(std::uint64_t circuit_seed) {
    // SplitMix64 gives W its own deterministic PRNG stream while the supplied
    // seed continues to select U exactly as it did before W was introduced.
    std::uint64_t value = circuit_seed + 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

void append_inverse_gates(Circuit& circuit, std::size_t first_gate) {
    if (first_gate > circuit.gates.size()) {
        throw std::invalid_argument("inverse range starts beyond the circuit");
    }
    const std::vector<Gate> forward_gates(
        circuit.gates.begin() + static_cast<std::ptrdiff_t>(first_gate),
        circuit.gates.end());
    circuit.gates.reserve(circuit.gates.size() + 2 * forward_gates.size());
    for (auto gate = forward_gates.rbegin();
         gate != forward_gates.rend();
         ++gate) {
        switch (gate->kind) {
            case GateKind::H:
                circuit.gates.push_back(Gate::h(gate->q0));
                break;
            case GateKind::S:
                // S^dagger = S^3. Keeping the existing gate vocabulary also
                // keeps the Pauli kernel's single canonical Clifford path.
                circuit.gates.push_back(Gate::s(gate->q0));
                circuit.gates.push_back(Gate::s(gate->q0));
                circuit.gates.push_back(Gate::s(gate->q0));
                break;
            case GateKind::CNOT:
                circuit.gates.push_back(Gate::cnot(gate->q0, gate->q1));
                break;
            case GateKind::RZ:
                circuit.gates.push_back(Gate::rz(gate->q0, -gate->parameter));
                break;
            case GateKind::Depolarizing:
                throw std::invalid_argument(
                    "a noisy channel has no unitary inverse for an identity "
                    "circuit");
        }
    }
}

double zero_state_pauli_expectation(const Pauli& observable) {
    return observable.x == 0 ? 1.0 : 0.0;
}

Pauli remap_pauli(
    const Pauli& pauli, const std::array<int, 64>& local_index) {
    Pauli result;
    for (int qubit = 0; qubit < 64; ++qubit) {
        if (local_index[qubit] < 0) {
            continue;
        }
        const std::uint64_t global_bit = std::uint64_t{1} << qubit;
        const std::uint64_t local_bit =
            std::uint64_t{1} << local_index[qubit];
        if ((pauli.x & global_bit) != 0) {
            result.x |= local_bit;
        }
        if ((pauli.z & global_bit) != 0) {
            result.z |= local_bit;
        }
    }
    return result;
}

Gate remap_gate(const Gate& gate, const std::array<int, 64>& local_index) {
    const int q0 = local_index[gate.q0];
    switch (gate.kind) {
        case GateKind::H:
            return Gate::h(q0);
        case GateKind::S:
            return Gate::s(q0);
        case GateKind::CNOT:
            return Gate::cnot(q0, local_index[gate.q1]);
        case GateKind::RZ:
            return Gate::rz(q0, gate.parameter);
        case GateKind::Depolarizing:
            throw std::invalid_argument(
                "lightcone statevector does not support noisy channels");
    }
    throw std::invalid_argument("unknown gate kind");
}

Circuit extract_gate_lightcone(const Circuit& circuit) {
    std::uint64_t active = circuit.observable.x | circuit.observable.z;
    std::vector<bool> relevant(circuit.gates.size(), false);
    for (std::size_t offset = 0; offset < circuit.gates.size(); ++offset) {
        const std::size_t index = circuit.gates.size() - 1 - offset;
        const Gate& gate = circuit.gates[index];
        const std::uint64_t q0_bit = std::uint64_t{1} << gate.q0;
        std::uint64_t gate_mask = q0_bit;
        if (gate.kind == GateKind::CNOT) {
            gate_mask |= std::uint64_t{1} << gate.q1;
        }
        if ((active & gate_mask) == 0) {
            continue;
        }
        relevant[index] = true;
        active |= gate_mask;
    }

    std::array<int, 64> local_index;
    local_index.fill(-1);
    int local_qubits = 0;
    for (int qubit = 0; qubit < circuit.qubits; ++qubit) {
        if ((active & (std::uint64_t{1} << qubit)) != 0) {
            local_index[qubit] = local_qubits++;
        }
    }

    Circuit lightcone;
    lightcone.family = circuit.family;
    lightcone.name = circuit.name + "_lightcone";
    lightcone.qubits = local_qubits;
    lightcone.observable = remap_pauli(circuit.observable, local_index);
    for (std::size_t index = 0; index < circuit.gates.size(); ++index) {
        if (relevant[index]) {
            lightcone.gates.push_back(remap_gate(circuit.gates[index], local_index));
        }
    }
    return lightcone;
}

}  // namespace

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
    circuit.family = depolarizing_probability > 0.0
                         ? "clifford_t_depol"
                         : "clifford_t";
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
    ReferenceQuery reference_query;
    reference_query.model = circuit.family;
    reference_query.qubits = qubits;
    reference_query.layers = layers;
    reference_query.circuit_generation_version = 1;
    reference_query.circuit_seed = seed;
    reference_query.t_density = t_density;
    reference_query.depolarizing_probability = depolarizing_probability;
    reference_query.observable = circuit.observable;
    attach_stored_reference(circuit, reference_query);

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
    ReferenceQuery reference_query;
    reference_query.model = circuit.family;
    reference_query.qubits = qubits;
    reference_query.layers = steps;
    reference_query.circuit_generation_version = 1;
    reference_query.dt = dt;
    reference_query.coupling = coupling;
    reference_query.transverse_field = transverse_field;
    reference_query.longitudinal_field = longitudinal_field;
    reference_query.observable = circuit.observable;
    attach_stored_reference(circuit, reference_query);

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

Circuit make_identity_echo_circuit(
    int qubits,
    int forward_layers,
    const std::string& base_model,
    double l1_cutoff) {
    Circuit circuit = make_prefixed_identity_echo_circuit(
        qubits,
        forward_layers,
        0,
        base_model,
        20260715,
        l1_cutoff);
    circuit.reference = zero_state_pauli_expectation(circuit.observable);
    circuit.reference_method = "analytic_identity";
    circuit.reference_id = "analytic_identity_" + circuit.name;
    return circuit;
}

Circuit make_prefixed_identity_echo_circuit(
    int qubits,
    int forward_layers,
    int prefix_layers,
    const std::string& base_model,
    std::uint64_t circuit_seed,
    double l1_cutoff) {
    if (prefix_layers < 0) {
        throw std::invalid_argument("W depth must be nonnegative");
    }
    if (forward_layers > std::numeric_limits<int>::max() - prefix_layers) {
        throw std::invalid_argument("combined circuit depth is too large");
    }

    Circuit circuit;
    std::size_t prefix_gate_count = 0;
    if (base_model == "clifford_t") {
        circuit = make_clifford_t_brickwork(
            qubits,
            forward_layers,
            0.70,
            circuit_seed,
            0.0,
            l1_cutoff);
        if (prefix_layers > 0) {
            Circuit prefix = make_clifford_t_brickwork(
                qubits,
                prefix_layers,
                0.70,
                prefix_circuit_seed(circuit_seed),
                0.0,
                l1_cutoff);
            prefix_gate_count = prefix.gates.size();
            prefix.gates.insert(
                prefix.gates.end(), circuit.gates.begin(), circuit.gates.end());
            circuit.gates = std::move(prefix.gates);
        }
        circuit.family = "clifford_t_identity";
        circuit.name = "ct_identity_n" + std::to_string(qubits) + "_w" +
                       std::to_string(prefix_layers) + "_u" +
                       std::to_string(forward_layers) + "_total_d" +
                       std::to_string(prefix_layers + 2 * forward_layers) +
                       "_seed" + std::to_string(circuit_seed);
    } else if (base_model == "ising") {
        circuit = make_nonintegrable_ising(
            qubits,
            forward_layers,
            0.12,
            1.0,
            0.91,
            0.37,
            l1_cutoff);
        if (prefix_layers > 0) {
            Circuit prefix = make_nonintegrable_ising(
                qubits, prefix_layers, 0.12, 1.0, 0.91, 0.37, l1_cutoff);
            prefix_gate_count = prefix.gates.size();
            prefix.gates.insert(
                prefix.gates.end(), circuit.gates.begin(), circuit.gates.end());
            circuit.gates = std::move(prefix.gates);
        }
        circuit.family = "ising_identity";
        circuit.name = "ising_identity_n" + std::to_string(qubits) + "_w" +
                       std::to_string(prefix_layers) + "_u" +
                       std::to_string(forward_layers) + "_total_s" +
                       std::to_string(prefix_layers + 2 * forward_layers);
    } else {
        throw std::invalid_argument(
            "identity circuit base model must be clifford_t or ising");
    }

    append_inverse_gates(circuit, prefix_gate_count);
    ReferenceQuery reference_query;
    reference_query.model = circuit.family;
    reference_query.qubits = qubits;
    reference_query.layers = forward_layers;
    reference_query.circuit_generation_version = 2;
    reference_query.circuit_seed =
        base_model == "clifford_t"
            ? std::optional<std::uint64_t>{circuit_seed}
            : std::nullopt;
    if (base_model == "clifford_t") {
        reference_query.t_density = 0.70;
        reference_query.depolarizing_probability = 0.0;
    } else {
        reference_query.dt = 0.12;
        reference_query.coupling = 1.0;
        reference_query.transverse_field = 0.91;
        reference_query.longitudinal_field = 0.37;
    }
    reference_query.prefix_depth = static_cast<double>(prefix_layers);
    reference_query.observable = circuit.observable;
    attach_stored_reference(circuit, reference_query);
    return circuit;
}

ExactIdentityReference compute_prefixed_identity_reference(
    int qubits,
    int prefix_layers,
    const std::string& base_model,
    std::uint64_t circuit_seed) {
    if (prefix_layers < 0) {
        throw std::invalid_argument("W depth must be nonnegative");
    }
    if (prefix_layers == 0) {
        const Circuit identity = make_identity_echo_circuit(
            qubits, 1, base_model, 0.0);
        return {
            zero_state_pauli_expectation(identity.observable),
            static_cast<int>(std::popcount(
                identity.observable.x | identity.observable.z)),
        };
    }

    if (base_model == "clifford_t") {
        const Circuit prefix = make_clifford_t_brickwork(
            qubits,
            prefix_layers,
            0.70,
            prefix_circuit_seed(circuit_seed),
            0.0,
            0.0);
        const Circuit lightcone = extract_gate_lightcone(prefix);
        return {
            statevector_expectation(lightcone),
            lightcone.qubits,
        };
    }
    if (base_model == "ising") {
        if (qubits < 2 || qubits > 64) {
            throw std::invalid_argument(
                "Ising circuits require between 2 and 64 qubits");
        }
        const int observable_start = std::clamp(qubits / 2 - 1, 0, qubits - 2);
        const int left = std::max(0, observable_start - prefix_layers);
        const int right = std::min(
            qubits - 1, observable_start + 1 + prefix_layers);
        const int lightcone_qubits = right - left + 1;
        Circuit lightcone = make_nonintegrable_ising(
            lightcone_qubits,
            prefix_layers,
            0.12,
            1.0,
            0.91,
            0.37,
            0.0);
        lightcone.observable.x = 0;
        lightcone.observable.z =
            (std::uint64_t{1} << (observable_start - left)) |
            (std::uint64_t{1} << (observable_start + 1 - left));
        return {
            statevector_expectation(lightcone),
            lightcone_qubits,
        };
    }
    throw std::invalid_argument(
        "identity circuit base model must be clifford_t or ising");
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
    if (model == "clifford_t_identity") {
        return make_identity_echo_circuit(
            qubits, layers, "clifford_t", l1_cutoff);
    }
    if (model == "ising_identity") {
        return make_identity_echo_circuit(
            qubits, layers, "ising", l1_cutoff);
    }
    throw std::invalid_argument(
        "model must be clifford_t, ising, clifford_t_depol, "
        "clifford_t_identity, or ising_identity");
}

std::vector<Circuit> benchmark_circuits() {
    return {
        make_configured_circuit(20, 12, "clifford_t", 0.00625),
        make_configured_circuit(20, 12, "ising", 0.00005),
        make_configured_circuit(20, 12, "clifford_t_depol", 0.034),
    };
}

}  // namespace pauli_bench
