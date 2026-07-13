#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace pb {
enum class GateKind { H, CNOT, RZ };
struct Gate { GateKind kind; int q0; int q1; double theta; };
struct Circuit {
 std::string name;
 int n;
 std::uint64_t observable_z_mask;
 std::vector<Gate> gates;
 int rz_count;
 double theta;
 double reference_value;
};
struct Pauli { std::uint64_t x{0}, z{0}; bool operator==(const Pauli&) const = default; };
struct PauliHash {
 std::size_t operator()(const Pauli& p) const noexcept {
  auto h=p.x*0x9E3779B185EBCA87ULL;
  h^=p.z+0x9E3779B97F4A7C15ULL+(h<<6)+(h>>2);
  return static_cast<std::size_t>(h);
 }
};
struct Result {
 std::string suite, case_name, method, status;
 int n{0}, depth{0}, rz_count{0};
 double theta{0}, seconds{0}, estimate{0}, exact_value{0}, abs_error{0}, std_error{0}, variance{0};
 std::uint64_t samples{0}, generated{0}, peak_terms{0}, terminal_paths{0};
};
}
