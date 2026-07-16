#pragma once

#include "pauli_bench/types.hpp"

#include <cstdint>
#include <vector>

namespace pauli_bench {

int local_label(const Pauli& p, int qubit) noexcept;
void set_local_label(Pauli& p, int qubit, int label) noexcept;

void apply_h_adjoint(Pauli& p, double& coefficient, int qubit) noexcept;
void apply_s_adjoint(Pauli& p, double& coefficient, int qubit) noexcept;
void apply_cnot_adjoint(Pauli& p, double& coefficient, int control, int target) noexcept;
void apply_depolarizing_adjoint(Pauli& p, double& coefficient, int qubit, double probability) noexcept;

class PauliKernel {
  public:
    PauliKernel();

    void apply_deterministic_in_place(Frontier& frontier, const Gate& gate) const;
    void apply_rz_and_merge(Frontier& frontier, const Gate& gate);

    [[nodiscard]] std::size_t merge_entries() const noexcept;
    [[nodiscard]] std::size_t merge_capacity() const noexcept;

  private:
    static constexpr double maximum_load_factor_ = 0.86;

    void begin_merge(std::size_t expected_entries);
    void ensure_capacity(std::size_t expected_entries);
    void grow(std::size_t new_capacity);
    void add(const Pauli& pauli, double coefficient);
    void insert_without_growth(const Pauli& pauli, double coefficient);

    std::vector<Pauli> keys_;
    std::vector<double> values_;
    std::vector<std::uint32_t> generations_;
    std::vector<std::uint32_t> occupied_;
    std::uint32_t generation_{1};
    std::size_t size_{};
    std::size_t mask_{};
};

}  // namespace pauli_bench
