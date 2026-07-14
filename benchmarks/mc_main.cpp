#include "pauli_bench/v1.hpp"
#include "pauli_bench/monte_carlo.hpp"
#include <iomanip>
#include <iostream>
using namespace pbv1;
int main(int argc,char**argv){double budget=argc>1?std::stod(argv[1]):30.;int idx=argc>2?std::stoi(argv[2]):-1;auto cs=benchmark_circuits();for(int i=0;i<(int)cs.size();++i){if(idx>=0&&i!=idx)continue;auto r=run_merged_subtree_monte_carlo(cs[i],budget,8.0,20260715);std::cout<<cs[i].name<<','<<std::setprecision(15)<<r.seconds<<','<<r.memory_mb<<','<<r.samples<<','<<r.estimate<<','<<r.std_error<<','<<r.error<<','<<r.peak_terms<<'\n';}}
