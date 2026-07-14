#include "pauli_bench/v1.hpp"
#include "pauli_bench/truncation.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace pbv1;
int main(int argc,char**argv){
 double bfs_budget=argc>1?std::stod(argv[1]):120.;
 std::ofstream o("results/benchmark_retained.csv");
 o<<"family,case,method,status,runtime_s,memory_mb,error,estimate,reference\n";
 for(auto&c:benchmark_circuits()){
  auto sv=run_statevector(c,120.); double ref=sv.estimate;
  auto bfs=run_bfs_l1_truncated(c,bfs_budget,1e-4);
  for(auto&r:{sv,bfs}){
   r.reference=ref; r.error=(r.status=="ok"?std::abs(r.estimate-ref):NAN);
   o<<r.family<<','<<r.case_name<<','<<r.method<<','<<r.status<<','<<std::setprecision(15)<<r.seconds<<','<<r.memory_mb<<','<<r.error<<','<<r.estimate<<','<<ref<<'\n';
   std::cout<<r.family<<" | "<<r.method<<" | "<<r.status<<" | t="<<r.seconds<<" s | mem="<<r.memory_mb<<" MB | err="<<r.error<<"\n";
  }
 }
}
