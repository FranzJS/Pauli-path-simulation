#include "pauli_bench/types.hpp"

#include <cmath>

namespace pauli_bench {
namespace {

std::uint64_t mix64(std::uint64_t value) noexcept {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

std::uint8_t checked_qubit(int q) {
    return static_cast<std::uint8_t>(q);
}

}  // namespace

Gate Gate::h(int q) {
    return {GateKind::H, checked_qubit(q), 0, 0.0, 1.0, 0.0};
}

Gate Gate::s(int q) {
    return {GateKind::S, checked_qubit(q), 0, 0.0, 1.0, 0.0};
}

Gate Gate::cnot(int control, int target) {
    return {
        GateKind::CNOT,
        checked_qubit(control),
        checked_qubit(target),
        0.0,
        1.0,
        0.0,
    };
}

Gate Gate::rz(int q, double theta) {
    return {
        GateKind::RZ,
        checked_qubit(q),
        0,
        theta,
        std::cos(theta),
        std::sin(theta),
    };
}

Gate Gate::depolarizing(int q, double probability) {
    return {
        GateKind::Depolarizing,
        checked_qubit(q),
        0,
        probability,
        1.0,
        0.0,
    };
}

std::size_t PauliHash::operator()(const Pauli& p) const noexcept {
    const auto hx = mix64(p.x);
    const auto hz = mix64(p.z ^ 0x6a09e667f3bcc909ULL);
    return static_cast<std::size_t>(
        hx ^ (hz + 0x9e3779b97f4a7c15ULL + (hx << 6U) + (hx >> 2U)));
}

}  // namespace pauli_bench
