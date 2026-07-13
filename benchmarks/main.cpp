#include "pauli_bench/circuits.hpp"
#include "pauli_bench/methods.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace pb;
static void emit(std::ofstream&o,const Result&r){o<<r.suite<<','<<r.case_name<<','<<r.method<<','<<r.status<<','<<r.n<<','<<r.depth<<','<<r.rz_count<<','<<std::setprecision(17)<<r.theta<<','<<r.seconds<<','<<r.estimate<<','<<r.exact_value<<','<<r.abs_error<<','<<r.std_error<<','<<r.variance<<','<<r.samples<<','<<r.generated<<','<<r.peak_terms<<','<<r.terminal_paths<<'\n';std::cout<<r.suite<<" "<<r.case_name<<" "<<r.method<<" "<<r.status<<" t="<<r.seconds<<" est="<<r.estimate<<"\n";}
int main(int argc,char**argv){double budget=argc>1?std::stod(argv[1]):25.;std::ofstream out("results/benchmark.csv");out<<"suite,case,method,status,n,depth,rz_count,theta,seconds,estimate,exact_value,abs_error,std_error,variance,samples,generated,peak_terms,terminal_paths\n";
 auto run=[&](std::string suite,const Circuit&c){Result truth=run_exact_sparse(c,budget);truth.suite=suite;truth.theta=c.theta;truth.exact_value=c.reference_value;truth.abs_error=std::abs(truth.estimate-c.reference_value);emit(out,truth);double exact=c.reference_value;
  if(c.n<=20){auto r=run_statevector(c,budget);r.suite=suite;r.theta=c.theta;r.exact_value=exact;r.abs_error=std::abs(r.estimate-exact);emit(out,r);}auto d=run_dfs(c,budget);d.suite=suite;d.theta=c.theta;d.exact_value=exact;d.abs_error=std::abs(d.estimate-exact);emit(out,d);auto m=run_monte_carlo(c,budget);m.suite=suite;m.theta=c.theta;m.exact_value=exact;m.abs_error=std::abs(m.estimate-exact);emit(out,m);};
 for(auto&c:width_suite())run("width",c);for(auto&c:branching_suite())run("branching",c);for(auto&c:variance_suite())run("variance",c);
}
