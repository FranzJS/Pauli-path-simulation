#include "pauli_bench/v1.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace pbv1;
int main(int argc,char**argv){
 double budget=argc>1?std::stod(argv[1]):25.;
 std::ofstream o("results/benchmark_v1.csv");
 o<<"family,case,method,status,runtime_s,memory_mb,error,estimate,reference\n";
 for(auto&c:benchmark_circuits()){
  auto sv=run_statevector(c,budget); double ref=sv.estimate;
  for(auto r:{sv,run_exact_sparse(c,budget),run_dfs(c,budget),run_monte_carlo(c,budget)}){
   r.reference=ref; r.error=(r.status=="ok"?std::abs(r.estimate-ref):NAN);
   o<<r.family<<','<<r.case_name<<','<<r.method<<','<<r.status<<','<<std::setprecision(12)<<r.seconds<<','<<r.memory_mb<<','<<r.error<<','<<r.estimate<<','<<ref<<'\n';
   std::cout<<r.family<<" | "<<r.method<<" | "<<r.status<<" | t="<<r.seconds<<" s | mem="<<r.memory_mb<<" MB | err="<<r.error<<"\n";
  }
 }
}
