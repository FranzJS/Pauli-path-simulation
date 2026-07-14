#pragma once
#include <cstdint>
#include <string>
#include <vector>
namespace pbv1 {
enum class GateKind { H, S, CNOT, RZ, DEPOL };
struct Gate { GateKind kind; int q0; int q1; double theta; };
struct Pauli { std::uint64_t x{0}, z{0}; bool operator==(const Pauli&) const = default; };
struct PauliHash { size_t operator()(const Pauli& p) const noexcept { auto h=p.x*0x9E3779B185EBCA87ULL; h^=p.z+0x9E3779B97F4A7C15ULL+(h<<6)+(h>>2); return static_cast<size_t>(h); } };
struct Circuit { std::string family,name; int n; Pauli observable; std::vector<Gate> gates; double cutoff; double reference; };
struct Result { std::string family,case_name,method,status; double seconds{0},memory_mb{0},estimate{0},reference{0},error{0},std_error{0},cutoff{0}; std::uint64_t peak_terms{0},samples{0}; };
Circuit make_clifford_t_brickwork(int,int,double,std::uint64_t,double noise_p=0);
Circuit make_nonintegrable_ising(int,int,double,double,double,double);
std::vector<Circuit> benchmark_circuits();
Result run_statevector(const Circuit&,double budget_s);
}
