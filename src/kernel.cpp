#include "pauli_bench/kernel.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace pauli_bench {
namespace {

std::size_t next_power_of_two(std::size_t value) {
    if (value <= 16) {
        return 16;
    }
    --value;
    for (std::size_t shift = 1; shift < sizeof(std::size_t) * 8; shift <<= 1U) {
        value |= value >> shift;
    }
    return value + 1;
}

}  // namespace

int local_label(const Pauli& p, int qubit) noexcept {
    const auto bit = std::uint64_t{1} << qubit;
    return ((p.x & bit) != 0 ? 1 : 0) | ((p.z & bit) != 0 ? 2 : 0);
}

void set_local_label(Pauli& p, int qubit, int label) noexcept {
    const auto bit = std::uint64_t{1} << qubit;
    p.x = (p.x & ~bit) | ((label & 1) != 0 ? bit : 0);
    p.z = (p.z & ~bit) | ((label & 2) != 0 ? bit : 0);
}

void apply_h_adjoint(Pauli& p, double& coefficient, int qubit) noexcept {
    const auto bit = std::uint64_t{1} << qubit;
    const bool x = (p.x & bit) != 0;
    const bool z = (p.z & bit) != 0;
    if (x && z) {
        coefficient = -coefficient;
    }
    if (x != z) {
        p.x ^= bit;
        p.z ^= bit;
    }
}

void apply_s_adjoint(Pauli& p, double& coefficient, int qubit) noexcept {
    const auto bit = std::uint64_t{1} << qubit;
    const bool x = (p.x & bit) != 0;
    const bool z = (p.z & bit) != 0;
    if (x && !z) {
        coefficient = -coefficient;
    }
    if (x) {
        p.z ^= bit;
    }
}

void apply_cnot_adjoint(
    Pauli& p,
    double& coefficient,
    int control,
    int target) noexcept {
    const auto control_bit = std::uint64_t{1} << control;
    const auto target_bit = std::uint64_t{1} << target;

    const bool x_control = (p.x & control_bit) != 0;
    const bool x_target = (p.x & target_bit) != 0;
    const bool z_control = (p.z & control_bit) != 0;
    const bool z_target = (p.z & target_bit) != 0;

    if (x_control && z_target && !(x_target ^ z_control)) {
        coefficient = -coefficient;
    }
    if (x_control) {
        p.x ^= target_bit;
    }
    if (z_target) {
        p.z ^= control_bit;
    }
}

void apply_depolarizing_adjoint(
    Pauli& p,
    double& coefficient,
    int qubit,
    double probability) noexcept {
    if (local_label(p, qubit) != 0) {
        coefficient *= 1.0 - probability;
    }
}

PauliKernel::PauliKernel() {
    grow(16);
}

void PauliKernel::begin_merge(std::size_t expected_entries) {
    ensure_capacity(expected_entries);
    ++generation_;
    if (generation_ == 0) {
        std::fill(generations_.begin(), generations_.end(), 0);
        generation_ = 1;
    }
    occupied_.clear();
    size_ = 0;
}

void PauliKernel::ensure_capacity(std::size_t expected_entries) {
    const auto required = static_cast<std::size_t>(
        std::ceil(static_cast<double>(std::max<std::size_t>(expected_entries, 1)) /
                  maximum_load_factor_));
    if (keys_.size() < required) {
        grow(next_power_of_two(required));
    }
}

void PauliKernel::grow(std::size_t new_capacity) {
    new_capacity = next_power_of_two(new_capacity);
    if (new_capacity > std::numeric_limits<std::uint32_t>::max()) {
        throw std::runtime_error("merge table exceeds 32-bit slot index capacity");
    }

    std::vector<Pauli> old_keys;
    std::vector<double> old_values;
    old_keys.reserve(size_);
    old_values.reserve(size_);
    if (!keys_.empty()) {
        for (const auto slot : occupied_) {
            old_keys.push_back(keys_[slot]);
            old_values.push_back(values_[slot]);
        }
    }

    keys_.assign(new_capacity, {});
    values_.assign(new_capacity, 0.0);
    generations_.assign(new_capacity, 0);
    occupied_.clear();
    occupied_.reserve(new_capacity / 2);
    generation_ = 1;
    size_ = 0;
    mask_ = new_capacity - 1;

    for (std::size_t index = 0; index < old_keys.size(); ++index) {
        insert_without_growth(old_keys[index], old_values[index]);
    }
}

void PauliKernel::insert_without_growth(const Pauli& pauli, double coefficient) {
    std::size_t slot = PauliHash{}(pauli) & mask_;
    while (generations_[slot] == generation_) {
        if (keys_[slot] == pauli) {
            values_[slot] += coefficient;
            return;
        }
        slot = (slot + 1) & mask_;
    }

    generations_[slot] = generation_;
    keys_[slot] = pauli;
    values_[slot] = coefficient;
    occupied_.push_back(static_cast<std::uint32_t>(slot));
    ++size_;
}

void PauliKernel::add(const Pauli& pauli, double coefficient) {
    if (coefficient == 0.0) {
        return;
    }
    if (size_ + 1 > static_cast<std::size_t>(keys_.size() * maximum_load_factor_)) {
        grow(keys_.size() * 2);
    }
    insert_without_growth(pauli, coefficient);
}

void PauliKernel::apply_deterministic_in_place(
    Frontier& frontier,
    const Gate& gate) const {
    assert(gate.kind != GateKind::RZ);

    switch (gate.kind) {
        case GateKind::H:
            for (auto& term : frontier) {
                apply_h_adjoint(term.pauli, term.coefficient, gate.q0);
            }
            break;
        case GateKind::S:
            for (auto& term : frontier) {
                apply_s_adjoint(term.pauli, term.coefficient, gate.q0);
            }
            break;
        case GateKind::CNOT:
            for (auto& term : frontier) {
                apply_cnot_adjoint(term.pauli, term.coefficient, gate.q0, gate.q1);
            }
            break;
        case GateKind::Depolarizing:
            for (auto& term : frontier) {
                apply_depolarizing_adjoint(
                    term.pauli,
                    term.coefficient,
                    gate.q0,
                    gate.parameter);
            }
            break;
        case GateKind::RZ:
            break;
    }
}

void PauliKernel::apply_rz_and_merge(Frontier& frontier, const Gate& gate) {
    assert(gate.kind == GateKind::RZ);

    begin_merge(frontier.size());

    const double cosine = gate.cos_parameter;
    const double sine = gate.sin_parameter;

    for (const auto& term : frontier) {
        const int label = local_label(term.pauli, gate.q0);
        if (label != 1 && label != 3) {
            add(term.pauli, term.coefficient);
            continue;
        }

        add(term.pauli, term.coefficient * cosine);

        Pauli rotated = term.pauli;
        if (label == 1) {
            set_local_label(rotated, gate.q0, 3);
            add(rotated, -term.coefficient * sine);
        } else {
            set_local_label(rotated, gate.q0, 1);
            add(rotated, term.coefficient * sine);
        }
    }

    frontier.clear();
    frontier.reserve(size_);
    for (const auto slot : occupied_) {
        const double coefficient = values_[slot];
        if (coefficient != 0.0) {
            frontier.push_back({keys_[slot], coefficient});
        }
    }
}

std::size_t PauliKernel::merge_entries() const noexcept {
    return size_;
}

std::size_t PauliKernel::merge_capacity() const noexcept {
    return keys_.size();
}

}  // namespace pauli_bench
