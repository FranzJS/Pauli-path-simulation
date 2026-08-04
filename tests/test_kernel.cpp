#include "pauli_bench/bfs.hpp"
#include "pauli_bench/circuits.hpp"
#include "pauli_bench/hybrid.hpp"
#include "pauli_bench/kernel.hpp"
#include "pauli_bench/reference.hpp"
#include "pauli_bench/truncation.hpp"
#include "pauli_bench/types.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
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

void test_configured_circuits() {
    using namespace pauli_bench;

    const auto clifford = make_configured_circuit(
        6, 3, "clifford_t", 0.0125);
    require(clifford.qubits == 6, "configured Clifford qubit count");
    require_close(clifford.l1_cutoff, 0.0125, 0.0, "configured L1 cutoff");
    require(
        clifford.observable == Pauli{0b1101, 0b0010},
        "configured Clifford observable must be centered");
    require(
        !std::isfinite(clifford.reference),
        "non-benchmark configuration must not reuse a fixed reference");

    const auto ising = make_configured_circuit(6, 2, "ising", 0.001);
    require(
        ising.observable == Pauli{0, 0b1100},
        "configured Ising observable must use the centered bond");

    const auto fixed = make_configured_circuit(
        20, 12, "clifford_t_depol", 0.034);
    require_close(
        fixed.reference,
        0.0663217307060528,
        0.0,
        "fixed configuration must retain its reference");

    bool rejected_model = false;
    try {
        static_cast<void>(make_configured_circuit(6, 3, "unknown", 0.01));
    } catch (const std::invalid_argument&) {
        rejected_model = true;
    }
    require(rejected_model, "unknown circuit model must be rejected");
}

void test_support_budget_truncation() {
    using namespace pauli_bench;

    Frontier frontier{
        {Pauli{1, 0}, 0.1},
        {Pauli{2, 0}, -0.2},
        {Pauli{3, 0}, 0.3},
        {Pauli{4, 0}, -0.4},
    };
    require(
        support_budget_size(frontier, 2, 0.2) == 2,
        "support budget K must combine threshold and cap");
    const auto removed = truncate_support_budget(frontier, 2, 0.2);
    require(removed == 2, "support truncation removed count");
    require(frontier.size() == 2, "support truncation must enforce cap");
    require(
        frontier[0].pauli.x == 3 && frontier[1].pauli.x == 4,
        "support truncation must retain the largest eligible terms");

    Frontier threshold_only{
        {Pauli{1, 0}, 0.1},
        {Pauli{2, 0}, 0.2},
        {Pauli{3, 0}, 0.3},
    };
    truncate_support_budget(threshold_only, 10, 0.2);
    require(
        threshold_only.size() == 2 && threshold_only[0].pauli.x == 2,
        "coefficient equal to threshold must be retained");

    const auto circuit = make_configured_circuit(
        6, 3, "clifford_t", 0.0);
    std::vector<std::size_t> schedule;
    const auto deterministic = run_bfs_support_truncated(
        circuit, 3, 0.0, 4, &schedule);
    require(
        deterministic.peak_post_truncation_terms <= 3,
        "deterministic support run must respect post-truncation cap");
    require(
        schedule.size() == deterministic.truncation_events,
        "support truncation must export one K per event");
    const auto automatic_pps = run_l1_optimal_pps_ht_support_budget(
        circuit, 1234, 3, 0.0, 4);
    require(
        automatic_pps.peak_post_truncation_terms <= 3,
        "automatic PPS must use support-budget K at every event");
    require(
        automatic_pps.truncation_events == deterministic.truncation_events,
        "automatic PPS and deterministic run must share event locations");

    KStrategyConfig scheduled_config;
    scheduled_config.strategy = KStrategy::Schedule;
    scheduled_config.schedule = schedule;
    const auto scheduled_deterministic = run_bfs_truncated(
        circuit, scheduled_config, 4);
    require_close(
        scheduled_deterministic.estimate,
        deterministic.estimate,
        0.0,
        "replaying an exported deterministic K schedule must reproduce estimate");

    const auto scheduled_pps = run_optimal_pps_ht(
        circuit, 1234, scheduled_config, 4);
    require(
        scheduled_pps.peak_post_truncation_terms <= 3,
        "scheduled PPS must consume the shared K strategy");

    KStrategyConfig invalid_pps_config;
    invalid_pps_config.strategy = KStrategy::L1Mass;
    invalid_pps_config.l1_cutoff = 0.01;
    bool rejected_l1_pps = false;
    try {
        static_cast<void>(run_optimal_pps_ht(
            circuit, 1234, invalid_pps_config, 4));
    } catch (const std::invalid_argument&) {
        rejected_l1_pps = true;
    }
    require(
        rejected_l1_pps,
        "PPS must reject deterministic-only L1 K determination");
}

void test_l1_optimal_pps_ht_truncation() {
    using namespace pauli_bench;

    for (std::uint64_t seed = 0; seed < 200; ++seed) {
        Frontier frontier{
            {Pauli{1, 0}, 1.0},
            {Pauli{2, 0}, 2.0},
            {Pauli{3, 0}, 100.0},
            {Pauli{4, 0}, 200.0},
        };
        std::mt19937_64 rng(seed);
        const auto stats = truncate_l1_optimal_pps_ht(frontier, 3, rng);

        require(frontier.size() == 3, "optimal PPS must retain exactly K terms");
        require(
            stats.heavy_terms == 2,
            "optimal PPS must derive its probability-one heavy hitters");
        require(
            stats.sampled_tail_terms == 1,
            "optimal PPS must use the remaining slot for randomized PPS");
        for (std::uint64_t key = 3; key <= 4; ++key) {
            const bool found = std::any_of(
                frontier.begin(),
                frontier.end(),
                [&](const Term& term) { return term.pauli.x == key; });
            require(found, "probability-one PPS coordinate must always be kept");
        }
    }

    Frontier zeros{
        {Pauli{1, 0}, 0.0},
        {Pauli{2, 0}, 0.0},
        {Pauli{3, 0}, 0.0},
        {Pauli{4, 0}, 0.0},
    };
    std::mt19937_64 zero_rng(7);
    const auto zero_stats = truncate_l1_optimal_pps_ht(zeros, 2, zero_rng);
    require(zeros.size() == 2, "zero-mass PPS must still honor budget K");
    require(
        zero_stats.max_importance_multiplier == 2.0,
        "zero-mass PPS must fall back to uniform inclusion probabilities");
    for (const auto& term : zeros) {
        require(
            std::isfinite(term.coefficient),
            "zero-mass PPS must not create non-finite coefficients");
    }

    // Exercise near-census systematic sampling with enough saturated entries
    // for accumulated probability rounding to cross an integer boundary.
    Frontier near_census;
    constexpr std::size_t heavy_terms = 23394;
    constexpr std::size_t population_terms = 27355;
    constexpr std::size_t retained_terms = 24000;
    near_census.reserve(population_terms);
    for (std::size_t index = 0; index < population_terms; ++index) {
        near_census.push_back({
            Pauli{static_cast<std::uint64_t>(index + 1), 0},
            index < population_terms - heavy_terms ? 1e-18 : 1.0});
    }
    std::mt19937_64 near_census_rng(81000000);
    const auto near_census_stats = truncate_l1_optimal_pps_ht(
        near_census, retained_terms, near_census_rng);
    require(
        near_census_stats.heavy_terms == heavy_terms,
        "near-census PPS must identify every saturated coordinate");
    require(
        static_cast<std::size_t>(std::count_if(
            near_census.begin(),
            near_census.end(),
            [](const Term& term) { return term.pauli.x > 3961; })) == heavy_terms,
        "near-census PPS must never omit saturated coordinates");
}

void test_l1_optimal_pps_ht_probabilities_and_unbiasedness() {
    using namespace pauli_bench;

    constexpr std::uint64_t trials = 12000;
    std::uint64_t selected_count[5]{};
    double coefficient_sum[5]{};
    for (std::uint64_t seed = 0; seed < trials; ++seed) {
        Frontier frontier{
            {Pauli{1, 0}, 1.0},
            {Pauli{2, 0}, 2.0},
            {Pauli{3, 0}, 3.0},
            {Pauli{4, 0}, 4.0},
        };
        std::mt19937_64 rng(seed);
        const auto stats = truncate_l1_optimal_pps_ht(frontier, 2, rng);

        require(frontier.size() == 2, "optimal PPS support must equal K");
        require(stats.heavy_terms == 0, "no coordinate should saturate here");
        require(stats.sampled_tail_terms == 2, "both slots must be randomized");
        double pass_sum = 0.0;
        for (const auto& term : frontier) {
            ++selected_count[term.pauli.x];
            coefficient_sum[term.pauli.x] += term.coefficient;
            pass_sum += term.coefficient;
            require_close(
                term.coefficient,
                5.0,
                1e-14,
                "PPS HT must equalize unsaturated adjusted magnitudes");
        }
        require_close(pass_sum, 10.0, 1e-14, "PPS HT total for positive inputs");
    }

    for (std::uint64_t key = 1; key <= 4; ++key) {
        const double expected_probability = static_cast<double>(key) / 5.0;
        const double frequency =
            static_cast<double>(selected_count[key]) /
            static_cast<double>(trials);
        require_close(
            frequency,
            expected_probability,
            0.02,
            "systematic PPS marginal inclusion probability");
        require_close(
            coefficient_sum[key] / static_cast<double>(trials),
            static_cast<double>(key),
            0.08,
            "optimal PPS HT coordinate mean must reproduce its input");
    }
}

}  // namespace

int main() {
    try {
        test_single_qubit_cliffords();
        test_cnot_generators();
        test_rz_branching();
        test_small_circuits_against_statevector();
        test_configured_circuits();
        test_support_budget_truncation();
        test_l1_optimal_pps_ht_truncation();
        test_l1_optimal_pps_ht_probabilities_and_unbiasedness();
    } catch (const std::exception& error) {
        std::cerr << "test failure: " << error.what() << '\n';
        return EXIT_FAILURE;
    }

    std::cout << "all tests passed\n";
    return EXIT_SUCCESS;
}
