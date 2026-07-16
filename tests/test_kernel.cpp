#include "pauli_bench/bfs.hpp"
#include "pauli_bench/kernel.hpp"
#include "pauli_bench/reference.hpp"
#include "pauli_bench/types.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void require_close(double actual, double expected, double tolerance, const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(
            message + ": expected " + std::to_string(expected) +
            ", got " + std::to_string(actual));
    }
}

void test_single_qubit_cliffords() {
    using namespace pauli_bench;

    {
        Pauli p{1, 0};
        double coefficient = 1.0;
        apply_h_adjoint(p, coefficient, 0);
        require(p == Pauli{0, 1}, "H: X must map to Z");
        require_close(coefficient, 1.0, 0.0, "H: X sign");
    }
    {
        Pauli p{0, 1};
        double coefficient = 1.0;
        apply_h_adjoint(p, coefficient, 0);
        require(p == Pauli{1, 0}, "H: Z must map to X");
    }
    {
        Pauli p{1, 1};
        double coefficient = 1.0;
        apply_h_adjoint(p, coefficient, 0);
        require(p == Pauli{1, 1}, "H: Y label");
        require_close(coefficient, -1.0, 0.0, "H: Y sign");
    }
    {
        Pauli p{1, 0};
        double coefficient = 1.0;
        apply_s_adjoint(p, coefficient, 0);
        require(p == Pauli{1, 1}, "S adjoint: X must map to Y label");
        require_close(coefficient, -1.0, 0.0, "S adjoint: X -> -Y");
    }
    {
        Pauli p{1, 1};
        double coefficient = 1.0;
        apply_s_adjoint(p, coefficient, 0);
        require(p == Pauli{1, 0}, "S adjoint: Y must map to X");
        require_close(coefficient, 1.0, 0.0, "S adjoint: Y sign");
    }
}

void test_cnot_generators() {
    using namespace pauli_bench;
    {
        Pauli p{1, 0};
        double coefficient = 1.0;
        apply_cnot_adjoint(p, coefficient, 0, 1);
        require(p.x == 3 && p.z == 0, "CNOT: Xc -> XcXt");
        require_close(coefficient, 1.0, 0.0, "CNOT Xc sign");
    }
    {
        Pauli p{2, 0};
        double coefficient = 1.0;
        apply_cnot_adjoint(p, coefficient, 0, 1);
        require(p.x == 2 && p.z == 0, "CNOT: Xt -> Xt");
    }
    {
        Pauli p{0, 1};
        double coefficient = 1.0;
        apply_cnot_adjoint(p, coefficient, 0, 1);
        require(p.x == 0 && p.z == 1, "CNOT: Zc -> Zc");
    }
    {
        Pauli p{0, 2};
        double coefficient = 1.0;
        apply_cnot_adjoint(p, coefficient, 0, 1);
        require(p.x == 0 && p.z == 3, "CNOT: Zt -> ZcZt");
    }
}

void test_rz_branching() {
    using namespace pauli_bench;
    Frontier frontier{{Pauli{1, 0}, 1.0}};
    PauliKernel kernel;
    kernel.apply_rz_and_merge(frontier, Gate::rz(0, std::numbers::pi / 4.0));
    require(frontier.size() == 2, "RZ X must branch to X and Y");

    double x = 0.0;
    double y = 0.0;
    for (const auto& term : frontier) {
        if (term.pauli == Pauli{1, 0}) {
            x = term.coefficient;
        } else if (term.pauli == Pauli{1, 1}) {
            y = term.coefficient;
        }
    }
    require_close(x, std::sqrt(0.5), 1e-14, "RZ X coefficient");
    require_close(y, -std::sqrt(0.5), 1e-14, "RZ Y coefficient");
}

void test_small_circuits_against_statevector() {
    using namespace pauli_bench;
    std::mt19937_64 rng(1234567);
    std::uniform_int_distribution<int> gate_choice(0, 3);
    std::uniform_int_distribution<int> qubit_choice(0, 4);

    for (int trial = 0; trial < 20; ++trial) {
        Circuit circuit;
        circuit.family = "test";
        circuit.name = "random_small";
        circuit.qubits = 5;
        circuit.observable.x = 1ULL << (trial % 5);
        circuit.observable.z = 1ULL << ((trial + 1) % 5);
        circuit.l1_cutoff = 0.0;
        circuit.reference_method = "statevector";

        for (int depth = 0; depth < 30; ++depth) {
            const int choice = gate_choice(rng);
            const int q0 = qubit_choice(rng);
            if (choice == 0) {
                circuit.gates.push_back(Gate::h(q0));
            } else if (choice == 1) {
                circuit.gates.push_back(Gate::s(q0));
            } else if (choice == 2) {
                int q1 = qubit_choice(rng);
                if (q1 == q0) {
                    q1 = (q1 + 1) % 5;
                }
                circuit.gates.push_back(Gate::cnot(q0, q1));
            } else {
                circuit.gates.push_back(Gate::rz(q0, 0.17 * (1 + (depth % 5))));
            }
        }

        circuit.reference = statevector_expectation(circuit);
        const auto result = run_bfs_l1_truncated(circuit, 4);
        require_close(
            result.estimate,
            circuit.reference,
            2e-11,
            "exact BFS must match statevector");
    }
}

}  // namespace

int main() {
    try {
        test_single_qubit_cliffords();
        test_cnot_generators();
        test_rz_branching();
        test_small_circuits_against_statevector();
    } catch (const std::exception& error) {
        std::cerr << "test failure: " << error.what() << '\n';
        return EXIT_FAILURE;
    }

    std::cout << "all tests passed\n";
    return EXIT_SUCCESS;
}
