#include "pauli_bench/terminal_aware_mc.hpp"
#include <chrono>
#include <cmath>
#include <random>
namespace pbv1 {
using Clock=std::chrono::steady_clock;
static double elapsed(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
static int label(const Pauli&p,int q){return (((p.x>>q)&1)?1:0)|(((p.z>>q)&1)?2:0);}
static void set_label(Pauli&p,int q,int l){auto m=1ULL<<q;p.x=(p.x&~m)|((l&1)?m:0);p.z=(p.z&~m)|((l&2)?m:0);}
static void apply_h(Pauli&p,double&w,int q){bool x=(p.x>>q)&1,z=(p.z>>q)&1;if(x&&z)w=-w;auto m=1ULL<<q;if(x!=z){p.x^=m;p.z^=m;}}
static void apply_s(Pauli&p,double&w,int q){bool x=(p.x>>q)&1,z=(p.z>>q)&1;if(x&&!z)w=-w;if(x)p.z^=1ULL<<q;}
static void apply_cnot(Pauli&p,double&w,int c,int t){bool xc=(p.x>>c)&1,xt=(p.x>>t)&1,zc=(p.z>>c)&1,zt=(p.z>>t)&1;if(xc&&zt&&(xt^zc^true))w=-w;if(xc)p.x^=1ULL<<t;if(zt)p.z^=1ULL<<c;}
static void deterministic(Pauli&p,double&w,const Gate&g){if(g.kind==GateKind::H)apply_h(p,w,g.q0);else if(g.kind==GateKind::S)apply_s(p,w,g.q0);else if(g.kind==GateKind::CNOT)apply_cnot(p,w,g.q0,g.q1);else if(g.kind==GateKind::DEPOL&&label(p,g.q0)!=0)w*=1.0-g.theta;}
static double proxy_score(const Circuit&c,Pauli p,int next_reverse_index,double rho,double floor){double a=1.0;const int total=static_cast<int>(c.gates.size());const int stop=std::min(total,next_reverse_index+64);for(int ri=next_reverse_index;ri<stop;++ri){const auto&g=c.gates[total-1-ri];if(g.kind==GateKind::RZ){int l=label(p,g.q0);if(l==1||l==3)a*=std::max(std::abs(std::cos(g.theta)),std::abs(std::sin(g.theta)));}else deterministic(p,a,g);}double terminal=std::pow(rho,__builtin_popcountll(p.x));return std::abs(a)*(floor+(1.0-floor)*terminal);}
Result run_terminal_aware_monte_carlo(const Circuit&c,double budget,std::uint64_t seed){Result r{c.family,c.name,"mc_terminal_aware","ok"};r.reference=c.reference;r.memory_mb=0.001;auto st=Clock::now();std::mt19937_64 rng(seed);std::uniform_real_distribution<double>u(0.0,1.0);double mean=0,m2=0;std::uint64_t n=0;const double rho=0.1;const double floor=c.family=="clifford_t_depol"?0.2:0.05;const int total=static_cast<int>(c.gates.size());while(elapsed(st)<budget){Pauli p=c.observable;double weight=1.0;for(int ri=0;ri<total;++ri){const auto&g=c.gates[total-1-ri];if(g.kind!=GateKind::RZ){deterministic(p,weight,g);continue;}int l=label(p,g.q0);if(l!=1&&l!=3)continue;Pauli p0=p,p1=p;double a0=std::cos(g.theta),a1=(l==1?-std::sin(g.theta):std::sin(g.theta));set_label(p1,g.q0,l==1?3:1);double s0=std::abs(a0)*proxy_score(c,p0,ri+1,rho,floor);double s1=std::abs(a1)*proxy_score(c,p1,ri+1,rho,floor);double den=s0+s1;if(!(den>0)){s0=std::abs(a0);s1=std::abs(a1);den=s0+s1;}double q0=s0/den;if(u(rng)<q0){p=p0;weight*=a0/q0;}else{p=p1;weight*=a1/(1.0-q0);}}
 double value=p.x==0?weight:0.0;++n;double d=value-mean;mean+=d/n;m2+=d*(value-mean);}r.seconds=elapsed(st);r.samples=n;r.estimate=mean;r.error=std::abs(mean-r.reference);if(n>1)r.std_error=std::sqrt((m2/(n-1))/n);return r;}
}
