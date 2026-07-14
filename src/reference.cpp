#include "pauli_bench/v1.hpp"
#include <chrono>
#include <cmath>
#include <complex>
#include <vector>
namespace pbv1 {
using C=std::complex<double>; using Clock=std::chrono::steady_clock;
static double elapsed(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
static void h(std::vector<C>&a,int q){auto s=1ULL<<q;double k=1/std::sqrt(2.);for(std::size_t b=0;b<a.size();b+=2*s)for(std::size_t i=0;i<s;++i){auto x=a[b+i],y=a[b+i+s];a[b+i]=(x+y)*k;a[b+i+s]=(x-y)*k;}}
static void sg(std::vector<C>&a,int q){auto m=1ULL<<q;for(std::size_t i=0;i<a.size();++i)if(i&m)a[i]*=C(0,1);}
static void rz(std::vector<C>&a,int q,double t){auto m=1ULL<<q;auto p0=std::exp(C(0,-t/2)),p1=std::exp(C(0,t/2));for(std::size_t i=0;i<a.size();++i)a[i]*=(i&m)?p1:p0;}
static void cx(std::vector<C>&a,int c,int t){auto cm=1ULL<<c,tm=1ULL<<t;for(std::size_t i=0;i<a.size();++i)if((i&cm)&&!(i&tm))std::swap(a[i],a[i|tm]);}
static C expectation(const std::vector<C>&a,const Pauli&p){C v=0;int yp=__builtin_popcountll(p.x&p.z)&3;C phase[4]={C(1,0),C(0,1),C(-1,0),C(0,-1)};for(std::size_t i=0;i<a.size();++i){auto j=i^p.x;double sign=__builtin_parityll(i&p.z)?-1:1;v+=std::conj(a[i])*a[j]*sign*phase[yp];}return v;}
Result run_statevector(const Circuit&c,double budget){Result r{c.family,c.name,"statevector","ok"};r.reference=c.reference;if(c.family=="clifford_t_depol"){r.status="not_applicable";return r;}auto st=Clock::now();std::vector<C>a(1ULL<<c.n);a[0]=1.;for(auto&g:c.gates){if(g.kind==GateKind::H)h(a,g.q0);else if(g.kind==GateKind::S)sg(a,g.q0);else if(g.kind==GateKind::RZ)rz(a,g.q0,g.theta);else if(g.kind==GateKind::CNOT)cx(a,g.q0,g.q1);if(elapsed(st)>budget){r.status="time_cap";break;}}if(r.status=="ok")r.estimate=expectation(a,c.observable).real();r.seconds=elapsed(st);r.memory_mb=a.size()*sizeof(C)/1048576.0;r.error=r.status=="ok"?std::abs(r.estimate-r.reference):NAN;return r;}
}
