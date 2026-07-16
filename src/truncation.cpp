#include "pauli_bench/truncation.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <unordered_map>
#include <vector>
namespace pbv1 {
using Clock=std::chrono::steady_clock; using Map=std::unordered_map<Pauli,double,PauliHash>;
static double elapsed(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
static int label(const Pauli&p,int q){return (((p.x>>q)&1)?1:0)|(((p.z>>q)&1)?2:0);} static void set_label(Pauli&p,int q,int l){auto m=1ULL<<q;p.x=(p.x&~m)|((l&1)?m:0);p.z=(p.z&~m)|((l&2)?m:0);}
static void ch(Pauli&p,double&w,int q){int l=label(p,q);if(l==1)set_label(p,q,2);else if(l==2)set_label(p,q,1);else if(l==3)w=-w;}
static void cs(Pauli&p,double&w,int q){int l=label(p,q);if(l==1){set_label(p,q,3);w=-w;}else if(l==3)set_label(p,q,1);}
struct CE{int a,b,s;};
static const CE& centry(int a,int b){static CE T[4][4];static bool init=false;if(!init){using C=std::complex<double>;C I[2][2]={{1,0},{0,1}},X[2][2]={{0,1},{1,0}},Z[2][2]={{1,0},{0,-1}},Y[2][2]={{0,C(0,-1)},{C(0,1),0}};C(*P[4])[2]={I,X,Z,Y};C U[4][4]={{1,0,0,0},{0,1,0,0},{0,0,0,1},{0,0,1,0}};for(int a0=0;a0<4;++a0)for(int b0=0;b0<4;++b0){C A[4][4]{};for(int i=0;i<2;++i)for(int j=0;j<2;++j)for(int k=0;k<2;++k)for(int l=0;l<2;++l)A[2*i+k][2*j+l]=P[a0][i][j]*P[b0][k][l];C B[4][4]{};for(int i=0;i<4;++i)for(int j=0;j<4;++j)for(int k=0;k<4;++k)for(int m=0;m<4;++m)B[i][j]+=std::conj(U[k][i])*A[k][m]*U[m][j];for(int ao=0;ao<4;++ao)for(int bo=0;bo<4;++bo){C tr=0.;for(int i=0;i<2;++i)for(int j=0;j<2;++j)for(int k=0;k<2;++k)for(int l=0;l<2;++l)tr+=std::conj(P[ao][i][j]*P[bo][k][l])*B[2*i+k][2*j+l];tr/=4.;if(std::abs(std::abs(tr.real())-1)<1e-9)T[a0][b0]={ao,bo,tr.real()>0?1:-1};}}init=true;}return T[a][b];}
static void ccx(Pauli&p,double&w,int c,int t){auto e=centry(label(p,c),label(p,t));set_label(p,c,e.a);set_label(p,t,e.b);w*=e.s;}
static void advance(const Map&cur,Map&nxt,const Gate&g){nxt.clear();nxt.reserve(cur.size()*2+1);for(auto const&[p0,w0]:cur){if(g.kind==GateKind::RZ){int l=label(p0,g.q0);if(l==1||l==3){double co=std::cos(g.theta),si=std::sin(g.theta);Pauli p2=p0;nxt[p0]+=w0*co;if(l==1){set_label(p2,g.q0,3);nxt[p2]+=w0*(-si);}else{set_label(p2,g.q0,1);nxt[p2]+=w0*si;}}else nxt[p0]+=w0;}else if(g.kind==GateKind::DEPOL){double w=w0;if(label(p0,g.q0)!=0)w*=1.0-g.theta;nxt[p0]+=w;}else{Pauli p=p0;double w=w0;if(g.kind==GateKind::H)ch(p,w,g.q0);else if(g.kind==GateKind::S)cs(p,w,g.q0);else ccx(p,w,g.q0,g.q1);nxt[p]+=w;}}}
static void truncate(Map&m,double frac){if(frac<=0||m.empty())return;std::vector<std::pair<Pauli,double>>v;v.reserve(m.size());double l1=0;for(auto&e:m){v.push_back(e);l1+=std::abs(e.second);}std::sort(v.begin(),v.end(),[](auto&a,auto&b){return std::abs(a.second)<std::abs(b.second);});double budget=frac*l1,used=0;size_t cut=0;while(cut<v.size()&&used+std::abs(v[cut].second)<=budget){used+=std::abs(v[cut].second);++cut;}if(!cut)return;m.clear();m.reserve((v.size()-cut)*2);for(size_t i=cut;i<v.size();++i)m.emplace(v[i]);}
static double finish(const Map&m){double v=0;for(auto const&[p,w]:m)if(p.x==0)v+=w;return v;}
Result run_bfs_l1_truncated(const Circuit&c,double budget,double cutoff,int interval){Result r{c.family,c.name,"bfs_l1","ok"};r.cutoff=cutoff;r.reference=c.reference;auto st=Clock::now();Map cur,nxt;cur[c.observable]=1.;std::uint64_t peak=1;int count=0;for(auto it=c.gates.rbegin();it!=c.gates.rend();++it){advance(cur,nxt,*it);cur.swap(nxt);peak=std::max<std::uint64_t>(peak,cur.size());if(it->kind==GateKind::RZ&&++count==interval){truncate(cur,cutoff);count=0;}if(elapsed(st)>budget){r.status="time_cap";break;}}if(count&&r.status=="ok")truncate(cur,cutoff);if(r.status=="ok")r.estimate=finish(cur);r.seconds=elapsed(st);r.peak_terms=peak;r.memory_mb=peak*48.0/1048576.0;r.error=r.status=="ok"?std::abs(r.estimate-r.reference):NAN;return r;}
}
