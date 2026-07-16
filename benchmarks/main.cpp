#include "pauli_bench/v1.hpp"
#include "pauli_bench/truncation.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace pbv1;
int main(int argc,char**argv){double budget=argc>1?std::stod(argv[1]):120.;std::ofstream o("results/benchmark_final.csv");o<<"family,case,method,status,runtime_s,memory_mb,cutoff,error,estimate,reference,peak_terms\n";for(auto&c:benchmark_circuits()){if(c.family!="clifford_t_depol"){auto sv=run_statevector(c,budget);o<<sv.family<<','<<sv.case_name<<','<<sv.method<<','<<sv.status<<','<<std::setprecision(15)<<sv.seconds<<','<<sv.memory_mb<<",,"<<sv.error<<','<<sv.estimate<<','<<sv.reference<<",\n";}auto r=run_bfs_l1_truncated(c,budget,c.cutoff,4);o<<r.family<<','<<r.case_name<<','<<r.method<<','<<r.status<<','<<r.seconds<<','<<r.memory_mb<<','<<r.cutoff<<','<<r.error<<','<<r.estimate<<','<<r.reference<<','<<r.peak_terms<<'\n';std::cout<<c.name<<" bfs t="<<r.seconds<<" mem="<<r.memory_mb<<" err="<<r.error<<"\n";}}
