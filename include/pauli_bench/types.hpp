#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace pauli_bench {

enum class GateKind : std::uint8_t { H, S, CNOT, RZ, Depolarizing };

struct Gate {
    GateKind kind{};
    std::uint8_t q0{};
    std::uint8_t q1{};
    double parameter{};
    double cos_parameter{1.0};
    double sin_parameter{0.0};

    static Gate h(int q);
    static Gate s(int q);
    static Gate cnot(int control, int target);
    static Gate rz(int q, double theta);
    static Gate depolarizing(int q, double probability);
};

struct Pauli {
    std::uint64_t x{};
    std::uint64_t z{};

    friend bool operator==(const Pauli&, const Pauli&) = default;
};

struct PauliHash {
    std::size_t operator()(const Pauli& p) const noexcept;
};

struct Term {
    Pauli pauli{};
    double coefficient{};
};

using Frontier = std::vector<Term>;

struct Circuit {
    std::string family;
    std::string name;
    int qubits{};
    Pauli observable{};
    std::vector<Gate> gates;
    double l1_cutoff{};
    double reference{};
    std::string reference_method;
    std::string reference_id;
};

struct BfsDiagnostics {
    double estimate{};
    double error{};
    double runtime_seconds{};
    std::uint64_t peak_support_terms{};
    std::uint64_t peak_pre_truncation_terms{};
    std::uint64_t peak_post_truncation_terms{};
    std::uint64_t truncation_events{};
    std::uint64_t peak_vector_capacity_terms{};
    std::uint64_t peak_merge_entries{};
    std::uint64_t peak_merge_capacity_slots{};
};

}  // namespace pauli_bench
