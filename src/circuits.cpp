#include "pauli_bench/circuits.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
namespace pb {
Circuit make_tensor_rx_circuit(int n, int active, double theta) {
 Circuit c{"tensor_rx_n"+std::to_string(n)+"_k"+std::to_string(active), n, 0, {}, 0, theta, std::pow(std::cos(theta), active)};
 for(int q=0;q<active;++q) c.observable_z_mask |= (1ULL<<q);
 for(int q=0;q<active;++q){
  c.gates.push_back({GateKind::H,q,0,0});
  c.gates.push_back({GateKind::RZ,q,0,theta});
  c.gates.push_back({GateKind::H,q,0,0});
  ++c.rz_count;
 }
 return c;
}
std::vector<Circuit> width_suite(){
 std::vector<Circuit> v;
 for(int n:{12,16,20,24,32,48,64}) v.push_back(make_tensor_rx_circuit(n,12,std::numbers::pi/8));
 return v;
}
std::vector<Circuit> branching_suite(){
 std::vector<Circuit> v;
 for(int k:{4,8,12,16,18,20,22}) v.push_back(make_tensor_rx_circuit(std::max(24,k),k,std::numbers::pi/4));
 return v;
}
std::vector<Circuit> variance_suite(){
 std::vector<Circuit> v;
 for(double t:{0.05,0.15,std::numbers::pi/8,std::numbers::pi/4}) v.push_back(make_tensor_rx_circuit(20,18,t));
 return v;
}
}
