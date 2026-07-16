#include "pauli_bench/reference.hpp"

#include <cmath>
#include <complex>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace pauli_bench {
namespace {

using Complex = std::complex<double>;

void apply_h(std::vector<Complex>& state, int qubit) {
    const std::size_t stride = std::size_t{1} << qubit;
    const double scale = 1.0 / std::sqrt(2.0);
    for (std::size_t block = 0; block < state.size(); block += 2 * stride) {
        for (std::size_t offset = 0; offset < stride; ++offset) {
            const Complex low = state[block + offset];
            const Complex high = state[block + offset + stride];
            state[block + offset] = (low + high) * scale;
            state[block + offset + stride] = (low - high) * scale;
        }
    }
}

void apply_s(std::vector<Complex>& state, int qubit) {
    const std::size_t bit = std::size_t{1} << qubit;
    for (std::size_t index = 0; index < state.size(); ++index) {
        if ((index & bit) != 0) {
            state[index] *= Complex{0.0, 1.0};
        }
    }
}

void apply_rz(std::vector<Complex>& state, int qubit, double theta) {
    const std::size_t bit = std::size_t{1} << qubit;
    const Complex phase_zero = std::exp(Complex{0.0, -theta / 2.0});
    const Complex phase_one = std::exp(Complex{0.0, theta / 2.0});
    for (std::size_t index = 0; index < state.size(); ++index) {
        state[index] *= (index & bit) != 0 ? phase_one : phase_zero;
    }
}

void apply_cnot(std::vector<Complex>& state, int control, int target) {
    const std::size_t control_bit = std::size_t{1} << control;
    const std::size_t target_bit = std::size_t{1} << target;
    for (std::size_t index = 0; index < state.size(); ++index) {
        if ((index & control_bit) != 0 && (index & target_bit) == 0) {
            std::swap(state[index], state[index | target_bit]);
        }
    }
}

Complex pauli_expectation(const std::vector<Complex>& state, const Pauli& pauli) {
    Complex result{0.0, 0.0};
    const int y_count_mod_four =
        __builtin_popcountll(pauli.x & pauli.z) & 3;
    const Complex phase[4] = {
        {1.0, 0.0},
        {0.0, 1.0},
        {-1.0, 0.0},
        {0.0, -1.0},
    };

    for (std::size_t index = 0; index < state.size(); ++index) {
        const std::size_t flipped = index ^ pauli.x;
        const double sign = __builtin_parityll(index & pauli.z) ? -1.0 : 1.0;
        result += std::conj(state[index]) * state[flipped] * sign *
                  phase[y_count_mod_four];
    }
    return result;
}

}  // namespace

double statevector_expectation(const Circuit& circuit) {
    if (circuit.qubits < 0 || circuit.qubits >= 63) {
        throw std::invalid_argument("unsupported statevector width");
    }

    std::vector<Complex> state(std::size_t{1} << circuit.qubits);
    state.front() = 1.0;

    for (const auto& gate : circuit.gates) {
        switch (gate.kind) {
            case GateKind::H:
                apply_h(state, gate.q0);
                break;
            case GateKind::S:
                apply_s(state, gate.q0);
                break;
            case GateKind::CNOT:
                apply_cnot(state, gate.q0, gate.q1);
                break;
            case GateKind::RZ:
                apply_rz(state, gate.q0, gate.parameter);
                break;
            case GateKind::Depolarizing:
                throw std::invalid_argument(
                    "statevector reference does not support noisy channels");
        }
    }

    return pauli_expectation(state, circuit.observable).real();
}

}  // namespace pauli_bench
