#include "pauli_bench/v1.hpp"
#include <chrono>
#include <cmath>
#include <complex>
#include <vector>
namespace pbv1 {
using Clock=std::chrono::steady_clock;
static double elapsed(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
static void apply_h(std::vector<std::complex<double>>&a,int q){auto s=1ULL<<q;double k=1/std::sqrt(2.);for(std::size_t b=0;b<a.size();b+=2*s)for(std::size_t i=0;i<s;++i){auto x=a[b+i],y=a[b+i+s];a[b+i]=(x+y)*k;a[b+i+s]=(x-y)*k;}}
static void apply_s(std::vector<std::complex<double>>&a,int q){auto m=1ULL<<q;for(std::size_t i=0;i<a.size();++i)if(i&m)a[i]*=std::complex<double>(0,1);}
static void apply_rz(std::vector<std::complex<double>>&a,int q,double t){auto m=1ULL<<q;auto p0=std::exp(std::complex<double>(0,-t/2)),p1=std::exp(std::complex<double>(0,t/2));for(std::size_t i=0;i<a.size();++i)a[i]*=(i&m)?p1:p0;}
static void apply_cnot(std::vector<std::complex<double>>&a,int c,int t){auto cm=1ULL<<c,tm=1ULL<<t;for(std::size_t i=0;i<a.size();++i)if((i&cm)&&!(i&tm))std::swap(a[i],a[i|tm]);}
Result run_statevector(const Circuit&c,double budget){
 Result r{c.family,c.name,"statevector","ok"};auto st=Clock::now();std::vector<std::complex<double>>a(1ULL<<c.n);a[0]=1.;
 for(auto&g:c.gates){if(g.kind==GateKind::H)apply_h(a,g.q0);else if(g.kind==GateKind::S)apply_s(a,g.q0);else if(g.kind==GateKind::RZ)apply_rz(a,g.q0,g.theta);else apply_cnot(a,g.q0,g.q1);if(elapsed(st)>budget){r.status="time_cap";break;}}
 if(r.status=="ok")for(std::size_t i=0;i<a.size();++i)r.estimate+=(__builtin_parityll(i&c.observable_z_mask)?-1:1)*std::norm(a[i]);
 r.seconds=elapsed(st);r.memory_mb=a.size()*sizeof(std::complex<double>)/1048576.0;return r;
}
}
