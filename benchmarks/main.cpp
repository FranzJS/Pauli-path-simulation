#include "pauli_bench/v1.hpp"
#include "pauli_bench/truncation.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
using namespace pbv1;
int main(int argc,char**argv){
 double budget=argc>1?std::stod(argv[1]):25.;
 std::ofstream o("results/benchmark_truncated.csv");
 o<<"family,case,method,status,runtime_s,memory_mb,error,estimate,reference\n";
 for(auto&c:benchmark_circuits()){
  auto sv=run_statevector(c,budget); double ref=sv.estimate;
  std::vector<Result> rs={sv,run_exact_sparse(c,budget),run_bfs_l1_truncated(c,budget),run_bfs_memory_capped(c,budget,2.),run_bfs_memory_capped(c,budget,4.),run_bfs_memory_capped(c,budget,8.),run_dfs(c,budget),run_monte_carlo(c,budget)};
  for(auto&r:rs){
   r.reference=ref; r.error=(r.status=="ok"?std::abs(r.estimate-ref):NAN);
   o<<r.family<<','<<r.case_name<<','<<r.method<<','<<r.status<<','<<std::setprecision(12)<<r.seconds<<','<<r.memory_mb<<','<<r.error<<','<<r.estimate<<','<<ref<<'\n';
   std::cout<<r.family<<" | "<<r.method<<" | "<<r.status<<" | t="<<r.seconds<<" s | mem="<<r.memory_mb<<" MB | err="<<r.error<<"\n";
  }
 }
}
