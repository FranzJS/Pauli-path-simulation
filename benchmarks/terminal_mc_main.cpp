#include "pauli_bench/terminal_aware_mc.hpp"
#include "pauli_bench/v1.hpp"
#include <iomanip>
#include <iostream>
using namespace pbv1;
int main(int argc,char**argv){double budget=argc>1?std::stod(argv[1]):10.;for(auto&c:benchmark_circuits()){auto r=run_terminal_aware_monte_carlo(c,budget,20260715);std::cout<<c.name<<','<<std::setprecision(15)<<r.seconds<<','<<r.samples<<','<<r.estimate<<','<<r.std_error<<','<<r.error<<','<<r.memory_mb<<'\n';}}
