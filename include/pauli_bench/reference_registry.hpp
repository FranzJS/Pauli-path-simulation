#pragma once

#include "pauli_bench/types.hpp"

#include <cstdint>
#include <optional>
#include <string>

namespace pauli_bench {

// Every field that identifies a generated circuit and measured observable.
// Optional values are represented by JSON null in the registry.
struct ReferenceQuery {
    std::string model;
    int qubits{};
    int layers{};
    int circuit_generation_version{1};
    std::optional<std::uint64_t> circuit_seed;
    std::optional<double> t_density;
    std::optional<double> depolarizing_probability;
    std::optional<double> dt;
    std::optional<double> coupling;
    std::optional<double> transverse_field;
    std::optional<double> longitudinal_field;
    std::optional<double> prefix_depth;
    Pauli observable{};
};

struct StoredReference {
    std::string id;
    double value{};
    std::string method;
    std::string precision;
    std::optional<double> uncertainty;
};

// PAULI_REFERENCE_REGISTRY overrides the repository registry path at runtime.
std::string reference_registry_path();

// Returns nullopt when the registry is valid but has no exact configuration
// match. Missing, malformed, or ambiguous registries raise runtime_error.
std::optional<StoredReference> find_stored_reference(
    const ReferenceQuery& query);

// Fill Circuit's reference fields, or mark them unavailable when no record
// matches. This is the normal entry point used by circuit construction.
void attach_stored_reference(Circuit& circuit, const ReferenceQuery& query);

}  // namespace pauli_bench
