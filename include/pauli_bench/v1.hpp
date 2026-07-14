#pragma once
#include <cstdint>
#include <string>
#include <vector>
namespace pbv1 {
enum class GateKind { H, S, CNOT, RZ, DEPOL };
struct Gate { GateKind kind; int q0; int q1; double theta; };
struct Circuit { std::string family,name; int n; std::uint64_t observable_z_mask; std::vector<Gate> gates; bool noisy{false}; double noise_p{0}; };
struct Pauli { std::uint64_t x{0},z{0}; bool operator==(const Pauli&) const = default; };
struct PauliHash { size_t operator()(const Pauli&p) const noexcept { auto h=p.x*0x9E3779B185EBCA87ULL; h^=p.z+0x9E3779B97F4A7C15ULL+(h<<6)+(h>>2); return (size_t)h; } };
struct Result { std::string family,case_name,method,status; double seconds{0},memory_mb{0},estimate{0},reference{0},error{0},std_error{0},cutoff{0}; };
Circuit make_clifford_t_brickwork(int n,int layers,double t_density,std::uint64_t seed,double noise_p=0);
Circuit make_nonintegrable_ising(int n,int steps,double dt,double J,double hx,double hz);
std::vector<Circuit> benchmark_circuits();
Result run_statevector(const Circuit&,double,int trajectories=1,std::uint64_t seed=42);
}
