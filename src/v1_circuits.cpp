#include "pauli_bench/v1.hpp"
#include <numbers>
#include <random>
namespace pbv1 {
Circuit make_clifford_t_brickwork(int n,int layers,double p,std::uint64_t seed){
 Circuit c{"clifford_t","ct_n"+std::to_string(n)+"_d"+std::to_string(layers),n,(0xFULL<<(n/2-2)),{}};
 std::mt19937_64 rng(seed); std::uniform_real_distribution<double> u(0,1); std::uniform_int_distribution<int> cliff(0,2);
 for(int l=0;l<layers;++l){
  for(int q=0;q<n;++q){ int r=cliff(rng); if(r==0)c.gates.push_back({GateKind::H,q,0,0}); else if(r==1)c.gates.push_back({GateKind::S,q,0,0}); else {c.gates.push_back({GateKind::H,q,0,0});c.gates.push_back({GateKind::S,q,0,0});} if(u(rng)<p)c.gates.push_back({GateKind::RZ,q,0,std::numbers::pi/4}); }
  int start=l&1; for(int q=start;q+1<n;q+=2)c.gates.push_back({GateKind::CNOT,q,q+1,0});
 }
 return c;
}
Circuit make_nonintegrable_ising(int n,int steps,double dt,double J,double hx,double hz){
 Circuit c{"ising","ising_n"+std::to_string(n)+"_s"+std::to_string(steps),n,1ULL<<(n/2),{}};
 for(int s=0;s<steps;++s){
  for(int q=0;q+1<n;++q){c.gates.push_back({GateKind::CNOT,q,q+1,0});c.gates.push_back({GateKind::RZ,q+1,0,2*J*dt});c.gates.push_back({GateKind::CNOT,q,q+1,0});}
  for(int q=0;q<n;++q){c.gates.push_back({GateKind::H,q,0,0});c.gates.push_back({GateKind::RZ,q,0,2*hx*dt});c.gates.push_back({GateKind::H,q,0,0});c.gates.push_back({GateKind::RZ,q,0,2*hz*dt});}
 }
 return c;
}
std::vector<Circuit> benchmark_circuits(){ return {make_clifford_t_brickwork(14,12,0.60,20260715), make_nonintegrable_ising(14,4,0.12,1.0,0.91,0.37)}; }
}
