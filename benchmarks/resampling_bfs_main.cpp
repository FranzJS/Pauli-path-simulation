#include "pauli_bench/resampling_bfs.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
using namespace pbv1;
static ResamplingMethod parse(std::string s){if(s=="ordinary")return ResamplingMethod::OrdinaryMultinomial;if(s=="residual")return ResamplingMethod::ResidualMultinomial;if(s=="dependent")return ResamplingMethod::ResidualDependentRounding;throw std::runtime_error("method must be ordinary, residual, or dependent");}
int main(int argc,char**argv){if(argc<5){std::cerr<<"usage: pauli_resampling_benchmark METHOD CIRCUIT_INDEX RETAINED_MB REPS [SEED]\n";return 2;}auto method=parse(argv[1]);int circuit_index=std::stoi(argv[2]);double mb=std::stod(argv[3]);int reps=std::stoi(argv[4]);std::uint64_t seed=argc>5?std::stoull(argv[5]):20260715;auto circuits=benchmark_circuits();if(circuit_index<0||circuit_index>=static_cast<int>(circuits.size()))return 2;const auto&c=circuits[circuit_index];double mean=0,m2=0,runtime=0,det=0;std::uint64_t peak=0,post=0,events=0;for(int r=0;r<reps;++r){auto x=run_resampling_bfs(c,mb,method,seed+1000003ULL*r);double delta=x.result.estimate-mean;mean+=delta/(r+1);m2+=delta*(x.result.estimate-mean);runtime+=x.result.seconds;det+=x.average_deterministic_fraction;peak=std::max(peak,x.result.peak_terms);post=std::max(post,x.peak_post_terms);events+=x.resampling_events;}double variance=reps>1?m2/(reps-1):std::numeric_limits<double>::quiet_NaN();double stderr=reps>1?std::sqrt(variance/reps):std::numeric_limits<double>::quiet_NaN();std::cout<<"method,case,retained_mb,reps,mean,variance,stderr,abs_error,mean_runtime_s,peak_pre_terms,peak_pre_mb,peak_post_terms,avg_deterministic_fraction,resampling_events_per_run\n"<<resampling_method_name(method)<<','<<c.name<<','<<mb<<','<<reps<<','<<std::setprecision(15)<<mean<<','<<variance<<','<<stderr<<','<<std::abs(mean-c.reference)<<','<<runtime/reps<<','<<peak<<','<<peak*48.0/1048576.0<<','<<post<<','<<det/reps<<','<<static_cast<double>(events)/reps<<'\n';}
