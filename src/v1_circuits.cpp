#include "pauli_bench/v1.hpp"
#include <numbers>
#include <random>
namespace pbv1 {
Circuit make_clifford_t_brickwork(int n,int layers,double p,std::uint64_t seed,double noise){
 Circuit c{noise>0?"clifford_t_depol":"clifford_t",(noise>0?"ct_depol_":"ct_")+std::string("n")+std::to_string(n)+"_d"+std::to_string(layers),n,(0xFULL<<(n/2-2)),{},noise>0,noise};
 std::mt19937_64 rng(seed); std::uniform_real_distribution<double> u(0,1); std::uniform_int_distribution<int> cliff(0,2);
 for(int l=0;l<layers;++l){
  for(int q=0;q<n;++q){int r=cliff(rng);if(r==0)c.gates.push_back({GateKind::H,q,0,0});else if(r==1)c.gates.push_back({GateKind::S,q,0,0});else{c.gates.push_back({GateKind::H,q,0,0});c.gates.push_back({GateKind::S,q,0,0});}if(u(rng)<p)c.gates.push_back({GateKind::RZ,q,0,std::numbers::pi/4});}
  for(int q=l&1;q+1<n;q+=2)c.gates.push_back({GateKind::CNOT,q,q+1,0});
  if(noise>0)for(int q=0;q<n;++q)c.gates.push_back({GateKind::DEPOL,q,0,noise});
 }
 return c;
}
Circuit make_nonintegrable_ising(int n,int steps,double dt,double J,double hx,double hz){
 Circuit c{"ising","ising_n"+std::to_string(n)+"_s"+std::to_string(steps),n,1ULL<<(n/2),{}};
 for(int s=0;s<steps;++s){for(int q=0;q+1<n;++q){c.gates.push_back({GateKind::CNOT,q,q+1,0});c.gates.push_back({GateKind::RZ,q+1,0,2*J*dt});c.gates.push_back({GateKind::CNOT,q,q+1,0});}for(int q=0;q<n;++q){c.gates.push_back({GateKind::H,q,0,0});c.gates.push_back({GateKind::RZ,q,0,2*hx*dt});c.gates.push_back({GateKind::H,q,0,0});c.gates.push_back({GateKind::RZ,q,0,2*hz*dt});}}
 return c;
}
std::vector<Circuit> benchmark_circuits(){return {make_clifford_t_brickwork(21,16,0.70,20260715,0),make_nonintegrable_ising(21,12,0.12,1.0,0.91,0.37),make_clifford_t_brickwork(21,16,0.70,20260715,0.05)};}
}
