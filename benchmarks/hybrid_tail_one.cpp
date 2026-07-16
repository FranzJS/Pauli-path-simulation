#include "pauli_bench/hybrid_tail.hpp"
#include "pauli_bench/truncation.hpp"
#include "pauli_bench/v1.hpp"
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
using namespace pbv1;
namespace {
void print_result(const Result&r,double ratio,std::uint64_t seed,std::uint64_t post=0,std::uint64_t heavy=0,std::uint64_t tail=0,double tail_mass=0,std::uint64_t events=0){
 std::cout<<r.family<<','<<r.case_name<<','<<r.method<<','<<ratio<<','<<seed<<','<<r.status<<','<<std::setprecision(15)<<r.seconds<<','<<r.memory_mb<<','<<r.peak_terms<<','<<post<<','<<heavy<<','<<tail<<','<<tail_mass<<','<<r.estimate<<','<<r.reference<<','<<r.error<<','<<events<<'\n';
}
}
int main(int argc,char**argv){
 try{
  if(argc<2)throw std::runtime_error("mode must be bfs or hybrid");
  const std::string mode=argv[1];auto circuits=benchmark_circuits();
  if(mode=="bfs"){
   if(argc<4)throw std::runtime_error("usage: pauli_hybrid_tail_one bfs CIRCUIT_INDEX BUDGET_S");
   const int circuit_index=std::stoi(argv[2]);const double budget=std::stod(argv[3]);const auto&c=circuits.at(circuit_index);auto r=run_bfs_l1_truncated(c,budget,c.cutoff,4);print_result(r,0,0);return r.status=="ok"?0:3;
  }
  if(mode=="hybrid"){
   if(argc<8)throw std::runtime_error("usage: pauli_hybrid_tail_one hybrid CIRCUIT_INDEX TAIL_RATIO BUDGET_S SEED SUPPORT_CAP WEIGHTING");
   const int circuit_index=std::stoi(argv[2]);const double ratio=std::stod(argv[3]);const double budget=std::stod(argv[4]);const std::uint64_t seed=std::stoull(argv[5]);const std::uint64_t cap=std::stoull(argv[6]);const std::string weighting_arg=argv[7];
   const auto weighting=weighting_arg=="ht"?TailWeighting::HorvitzThompson:TailWeighting::OriginalCoefficient;const auto&c=circuits.at(circuit_index);auto d=run_hybrid_l1_tail(c,budget,c.cutoff,ratio,seed,4,weighting,cap);print_result(d.result,ratio,seed,d.peak_post_terms,d.max_heavy_terms,d.max_sampled_tail_terms,d.average_retained_tail_l1_fraction,d.truncation_events);return d.result.status=="ok"?0:3;
  }
  throw std::runtime_error("mode must be bfs or hybrid");
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
}
