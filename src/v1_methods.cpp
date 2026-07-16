#include "pauli_bench/v1.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <random>
#include <unordered_map>
#include <vector>
namespace pbv1 {
using Clock=std::chrono::steady_clock;
static double elapsed(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
static int label(const Pauli&p,int q){return (((p.x>>q)&1)?1:0)|(((p.z>>q)&1)?2:0);}
static void set_label(Pauli&p,int q,int l){auto m=1ULL<<q;p.x=(p.x&~m)|((l&1)?m:0);p.z=(p.z&~m)|((l&2)?m:0);}
static void conj_h(Pauli&p,double&w,int q){int l=label(p,q);if(l==1)set_label(p,q,2);else if(l==2)set_label(p,q,1);else if(l==3)w=-w;}
static void conj_s(Pauli&p,double&w,int q){int l=label(p,q);if(l==1){set_label(p,q,3);w=-w;}else if(l==3)set_label(p,q,1);}
struct CEntry{int a,b,s;};
static const CEntry& centry(int a,int b){
 static CEntry T[4][4]; static bool init=false;
 if(!init){using C=std::complex<double>;C I[2][2]={{1,0},{0,1}},X[2][2]={{0,1},{1,0}},Z[2][2]={{1,0},{0,-1}},Y[2][2]={{0,C(0,-1)},{C(0,1),0}};C(*P[4])[2]={I,X,Z,Y};C U[4][4]={{1,0,0,0},{0,1,0,0},{0,0,0,1},{0,0,1,0}};
  for(int a0=0;a0<4;++a0)for(int b0=0;b0<4;++b0){C A[4][4]{};for(int i=0;i<2;++i)for(int j=0;j<2;++j)for(int k=0;k<2;++k)for(int l=0;l<2;++l)A[2*i+k][2*j+l]=P[a0][i][j]*P[b0][k][l];C B[4][4]{};for(int i=0;i<4;++i)for(int j=0;j<4;++j)for(int k=0;k<4;++k)for(int m=0;m<4;++m)B[i][j]+=std::conj(U[k][i])*A[k][m]*U[m][j];for(int ao=0;ao<4;++ao)for(int bo=0;bo<4;++bo){C tr=0.;for(int i=0;i<2;++i)for(int j=0;j<2;++j)for(int k=0;k<2;++k)for(int l=0;l<2;++l)tr+=std::conj(P[ao][i][j]*P[bo][k][l])*B[2*i+k][2*j+l];tr/=4.;if(std::abs(std::abs(tr.real())-1)<1e-9)T[a0][b0]={ao,bo,tr.real()>0?1:-1};}}
  init=true;}return T[a][b];}
static void conj_cnot(Pauli&p,double&w,int c,int t){auto e=centry(label(p,c),label(p,t));set_label(p,c,e.a);set_label(p,t,e.b);w*=e.s;}
static bool z_only(const Pauli&p){return p.x==0;}
static void cliff(Pauli&p,double&w,const Gate&g){if(g.kind==GateKind::H)conj_h(p,w,g.q0);else if(g.kind==GateKind::S)conj_s(p,w,g.q0);else conj_cnot(p,w,g.q0,g.q1);}

Result run_exact_sparse(const Circuit&c,double budget,std::uint64_t cap){Result r{c.family,c.name,"exact_sparse","ok"};auto st=Clock::now();std::unordered_map<Pauli,double,PauliHash> cur,nxt;Pauli o;o.z=c.observable_z_mask;cur[o]=1.;std::uint64_t peak=1;
 for(auto it=c.gates.rbegin();it!=c.gates.rend();++it){nxt.clear();nxt.reserve(cur.size()*2+1);for(auto const&[p0,w0]:cur){if(it->kind==GateKind::RZ){int l=label(p0,it->q0);if(l==1||l==3){double co=std::cos(it->theta),si=std::sin(it->theta);Pauli p2=p0;nxt[p0]+=w0*co;if(l==1){set_label(p2,it->q0,3);nxt[p2]+=w0*(-si);}else{set_label(p2,it->q0,1);nxt[p2]+=w0*si;}}else nxt[p0]+=w0;}else{Pauli p=p0;double w=w0;cliff(p,w,*it);nxt[p]+=w;}}
  cur.swap(nxt);peak=std::max<std::uint64_t>(peak,cur.size());if(cur.size()>cap){r.status="memory_cap";break;}if(elapsed(st)>budget){r.status="time_cap";break;}}
 if(r.status=="ok")for(auto const&[p,w]:cur)if(z_only(p))r.estimate+=w;r.seconds=elapsed(st);r.memory_mb=peak*48.0/1048576.0;return r;}

struct DS{const Circuit*c;Clock::time_point st;double budget;std::uint64_t cap,paths=0,nodes=0;bool stop=false;double sum=0;};
static void rec(DS&s,int i,Pauli p,double w){if(s.stop)return;if((++s.nodes&0x3ffff)==0&&elapsed(s.st)>s.budget){s.stop=true;return;}if(i<0){if(z_only(p))s.sum+=w;if(++s.paths>=s.cap)s.stop=true;return;}auto&g=s.c->gates[i];if(g.kind==GateKind::RZ){int l=label(p,g.q0);if(l==1||l==3){double co=std::cos(g.theta),si=std::sin(g.theta);rec(s,i-1,p,w*co);Pauli p2=p;if(l==1){set_label(p2,g.q0,3);rec(s,i-1,p2,w*(-si));}else{set_label(p2,g.q0,1);rec(s,i-1,p2,w*si);}return;}}if(g.kind!=GateKind::RZ)cliff(p,w,g);rec(s,i-1,p,w);}
Result run_dfs(const Circuit&c,double budget,std::uint64_t cap){Result r{c.family,c.name,"dfs","ok"};DS s{&c,Clock::now(),budget,cap};Pauli p;p.z=c.observable_z_mask;rec(s,(int)c.gates.size()-1,p,1.);r.seconds=elapsed(s.st);r.estimate=s.sum;r.memory_mb=(c.gates.size()*sizeof(Pauli)*2.0)/1048576.0;if(s.stop)r.status=s.paths>=cap?"path_cap":"time_cap";return r;}

Result run_monte_carlo(const Circuit&c,double budget,std::uint64_t maxs,std::uint64_t seed){Result r{c.family,c.name,"monte_carlo","ok"};auto st=Clock::now();std::mt19937_64 rng(seed);std::uniform_real_distribution<double>u(0,1);long double sum=0;std::uint64_t S=0;
 while(S<maxs&&elapsed(st)<budget){std::uint64_t take=std::min<std::uint64_t>(4096,maxs-S);for(std::uint64_t k=0;k<take;++k){Pauli p;p.z=c.observable_z_mask;double w=1.;for(auto it=c.gates.rbegin();it!=c.gates.rend();++it){if(it->kind==GateKind::RZ){int l=label(p,it->q0);if(l==1||l==3){double a=std::abs(std::cos(it->theta)),b=std::abs(std::sin(it->theta)),z=a+b,pr=a/z;if(u(rng)<pr)w*=std::cos(it->theta)/pr;else{w*=(l==1?-std::sin(it->theta):std::sin(it->theta))/(1-pr);set_label(p,it->q0,l==1?3:1);}}}else cliff(p,w,*it);}if(z_only(p))sum+=w;}S+=take;}
 r.estimate=S?double(sum/S):0;r.seconds=elapsed(st);r.memory_mb=0.01;if(S<maxs)r.status="time_cap";return r;}

static void ah(std::vector<std::complex<double>>&a,int q){auto s=1ULL<<q;double k=1/std::sqrt(2.);for(size_t b=0;b<a.size();b+=2*s)for(size_t i=0;i<s;++i){auto x=a[b+i],y=a[b+i+s];a[b+i]=(x+y)*k;a[b+i+s]=(x-y)*k;}}
static void as(std::vector<std::complex<double>>&a,int q){auto m=1ULL<<q;for(size_t i=0;i<a.size();++i)if(i&m)a[i]*=std::complex<double>(0,1);}
static void arz(std::vector<std::complex<double>>&a,int q,double t){auto m=1ULL<<q;auto p0=std::exp(std::complex<double>(0,-t/2)),p1=std::exp(std::complex<double>(0,t/2));for(size_t i=0;i<a.size();++i)a[i]*=(i&m)?p1:p0;}
static void acx(std::vector<std::complex<double>>&a,int c,int t){auto cm=1ULL<<c,tm=1ULL<<t;for(size_t i=0;i<a.size();++i)if((i&cm)&&!(i&tm))std::swap(a[i],a[i|tm]);}
Result run_statevector(const Circuit&c,double budget){Result r{c.family,c.name,"statevector","ok"};auto st=Clock::now();std::vector<std::complex<double>>a(1ULL<<c.n);a[0]=1.;for(auto&g:c.gates){if(g.kind==GateKind::H)ah(a,g.q0);else if(g.kind==GateKind::S)as(a,g.q0);else if(g.kind==GateKind::RZ)arz(a,g.q0,g.theta);else acx(a,g.q0,g.q1);if(elapsed(st)>budget){r.status="time_cap";break;}}if(r.status=="ok"){for(size_t i=0;i<a.size();++i)r.estimate+=(__builtin_parityll(i&c.observable_z_mask)?-1:1)*std::norm(a[i]);}r.seconds=elapsed(st);r.memory_mb=a.size()*sizeof(std::complex<double>)/1048576.0;return r;}
}
