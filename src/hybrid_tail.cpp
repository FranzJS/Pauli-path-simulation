#include "pauli_bench/hybrid_tail.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <limits>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>
namespace pbv1 {
namespace {
using Map=std::unordered_map<Pauli,double,PauliHash>;
using Clock=std::chrono::steady_clock;
int label(const Pauli&p,int q){return (((p.x>>q)&1)?1:0)|(((p.z>>q)&1)?2:0);}
void set_label(Pauli&p,int q,int l){auto m=1ULL<<q;p.x=(p.x&~m)|((l&1)?m:0);p.z=(p.z&~m)|((l&2)?m:0);}
void apply_h(Pauli&p,double&w,int q){int l=label(p,q);if(l==1)set_label(p,q,2);else if(l==2)set_label(p,q,1);else if(l==3)w=-w;}
void apply_s(Pauli&p,double&w,int q){int l=label(p,q);if(l==1){set_label(p,q,3);w=-w;}else if(l==3)set_label(p,q,1);}
struct CnotEntry{int a=0,b=0,s=1;};
const CnotEntry& cnot_entry(int a,int b){
 static CnotEntry table[4][4]; static bool init=false;
 if(!init){using C=std::complex<double>;C I[2][2]={{1,0},{0,1}},X[2][2]={{0,1},{1,0}},Z[2][2]={{1,0},{0,-1}},Y[2][2]={{0,C(0,-1)},{C(0,1),0}};C(*P[4])[2]={I,X,Z,Y};C U[4][4]={{1,0,0,0},{0,1,0,0},{0,0,0,1},{0,0,1,0}};for(int a0=0;a0<4;++a0)for(int b0=0;b0<4;++b0){C A[4][4]{};for(int i=0;i<2;++i)for(int j=0;j<2;++j)for(int k=0;k<2;++k)for(int l=0;l<2;++l)A[2*i+k][2*j+l]=P[a0][i][j]*P[b0][k][l];C B[4][4]{};for(int i=0;i<4;++i)for(int j=0;j<4;++j)for(int k=0;k<4;++k)for(int m=0;m<4;++m)B[i][j]+=std::conj(U[k][i])*A[k][m]*U[m][j];for(int ao=0;ao<4;++ao)for(int bo=0;bo<4;++bo){C tr=0.;for(int i=0;i<2;++i)for(int j=0;j<2;++j)for(int k=0;k<2;++k)for(int l=0;l<2;++l)tr+=std::conj(P[ao][i][j]*P[bo][k][l])*B[2*i+k][2*j+l];tr/=4.;if(std::abs(std::abs(tr.real())-1)<1e-9)table[a0][b0]={ao,bo,tr.real()>0?1:-1};}}init=true;}
 return table[a][b];
}
void apply_cnot(Pauli&p,double&w,int c,int t){auto e=cnot_entry(label(p,c),label(p,t));set_label(p,c,e.a);set_label(p,t,e.b);w*=e.s;}
enum class AdvanceStatus { Ok, TimeCap, SupportCap };
AdvanceStatus advance(const Map&cur,Map&next,const Gate&g,Clock::time_point start,double budget_s,std::uint64_t support_cap){
 next.clear();
 const auto reserve_target=std::min<std::uint64_t>(support_cap,(g.kind==GateKind::RZ?2ULL:1ULL)*cur.size()+1ULL);
 next.reserve(static_cast<std::size_t>(reserve_target));
 std::uint64_t processed=0;
 for(const auto&[p0,w0]:cur){
  if(g.kind==GateKind::RZ){int l=label(p0,g.q0);if(l==1||l==3){double co=std::cos(g.theta),si=std::sin(g.theta);Pauli p2=p0;next[p0]+=w0*co;if(l==1){set_label(p2,g.q0,3);next[p2]+=w0*(-si);}else{set_label(p2,g.q0,1);next[p2]+=w0*si;}}else next[p0]+=w0;}
  else if(g.kind==GateKind::DEPOL){double w=w0;if(label(p0,g.q0)!=0)w*=1.0-g.theta;next[p0]+=w;}
  else{Pauli p=p0;double w=w0;if(g.kind==GateKind::H)apply_h(p,w,g.q0);else if(g.kind==GateKind::S)apply_s(p,w,g.q0);else apply_cnot(p,w,g.q0,g.q1);next[p]+=w;}
  if(next.size()>support_cap)return AdvanceStatus::SupportCap;
  if((++processed&2047ULL)==0&&std::chrono::duration<double>(Clock::now()-start).count()>budget_s)return AdvanceStatus::TimeCap;
 }
 return AdvanceStatus::Ok;
}
std::vector<double> inclusion_probabilities(const std::vector<double>&w,std::size_t target){
 const std::size_t n=w.size();std::vector<double>pi(n,0.0);if(target==0||n==0)return pi;if(target>=n){std::fill(pi.begin(),pi.end(),1.0);return pi;}
 std::vector<long double>prefix(n+1,0.0L);for(std::size_t i=0;i<n;++i)prefix[i+1]=prefix[i]+static_cast<long double>(w[i]);
 std::size_t certain=0;
 while(certain<target){
  const std::size_t remaining_count=n-certain;const long double remaining=prefix[remaining_count];
  if(!(remaining>0))throw std::runtime_error("nonpositive PPS tail mass");
  const long double lambda=static_cast<long double>(target-certain)/remaining;
  const std::size_t largest_remaining=remaining_count-1;
  if(lambda*static_cast<long double>(w[largest_remaining])<1.0L-1e-15L)break;
  pi[largest_remaining]=1.0;++certain;
 }
 const std::size_t fractional_count=n-certain;
 if(certain<target){const long double lambda=static_cast<long double>(target-certain)/prefix[fractional_count];for(std::size_t i=0;i<fractional_count;++i)pi[i]=static_cast<double>(std::min<long double>(1.0L,lambda*static_cast<long double>(w[i])));}
 double sum=0;for(double value:pi)sum+=value;double diff=static_cast<double>(target)-sum;
 if(std::abs(diff)>1e-12){for(std::size_t i=0;i<fractional_count&&std::abs(diff)>1e-12;++i){double room=diff>0?1.0-pi[i]:pi[i];double delta=std::copysign(std::min(std::abs(diff),room),diff);pi[i]+=delta;diff-=delta;}}
 return pi;
}
std::vector<std::uint8_t> systematic_sample(const std::vector<double>&pi,std::mt19937_64&rng,std::size_t target){
 std::vector<std::uint8_t>selected(pi.size(),0);if(target==0)return selected;
 std::uniform_real_distribution<double>uniform(0.0,1.0);long double threshold=uniform(rng),cumulative=0.0L;std::size_t count=0;
 for(std::size_t i=0;i<pi.size()&&count<target;++i){cumulative+=static_cast<long double>(pi[i]);if(threshold<cumulative+1e-15L){selected[i]=1;++count;threshold+=1.0L;}}
 if(count!=target){for(std::size_t i=pi.size();i-->0&&count<target;)if(!selected[i]&&pi[i]>0){selected[i]=1;++count;}}
 if(count!=target)throw std::runtime_error("systematic PPS population mismatch");
 return selected;
}
struct CompressionStats{std::size_t heavy=0,sampled=0;double sampled_tail_l1_fraction=0;};
CompressionStats compress(Map&m,double cutoff,double tail_ratio,std::mt19937_64&rng,TailWeighting weighting){
 struct Entry{Pauli p;double w;double a;};std::vector<Entry>v;v.reserve(m.size());double l1=0;
 for(const auto&[p,w]:m)if(w!=0){double a=std::abs(w);v.push_back({p,w,a});l1+=a;}
 std::sort(v.begin(),v.end(),[](const Entry&a,const Entry&b){if(a.a!=b.a)return a.a<b.a;if(a.p.x!=b.p.x)return a.p.x<b.p.x;return a.p.z<b.p.z;});
 double budget=cutoff*l1,used=0;std::size_t cut=0;while(cut<v.size()&&used+v[cut].a<=budget){used+=v[cut].a;++cut;}
 const std::size_t heavy=v.size()-cut;std::size_t requested=heavy?static_cast<std::size_t>(std::ceil(tail_ratio*static_cast<double>(heavy))):0;
 std::vector<std::size_t>positive_tail;positive_tail.reserve(cut);std::vector<double>weights;weights.reserve(cut);for(std::size_t i=0;i<cut;++i)if(v[i].a>0){positive_tail.push_back(i);weights.push_back(v[i].a);}
 const std::size_t sampled=std::min(requested,positive_tail.size());Map out;out.reserve((heavy+sampled)*2+1);for(std::size_t i=cut;i<v.size();++i)out.emplace(v[i].p,v[i].w);
 double sampled_mass=0;
 if(sampled){auto pi=inclusion_probabilities(weights,sampled);auto selected=systematic_sample(pi,rng,sampled);for(std::size_t j=0;j<selected.size();++j)if(selected[j]){auto&entry=v[positive_tail[j]];out.emplace(entry.p,weighting==TailWeighting::HorvitzThompson?entry.w/pi[j]:entry.w);sampled_mass+=entry.a;}}
 m.swap(out);return {heavy,sampled,used>0?sampled_mass/used:0.0};
}
double finish(const Map&m){double value=0;for(const auto&[p,w]:m)if(p.x==0)value+=w;return value;}
}
HybridTailDiagnostics run_hybrid_l1_tail(const Circuit&c,double budget_s,double cutoff,double tail_ratio,std::uint64_t seed,int interval,TailWeighting weighting,std::uint64_t support_cap){
 HybridTailDiagnostics d;d.result={c.family,c.name,(weighting==TailWeighting::HorvitzThompson?"hybrid_l1_tail_ht_":"hybrid_l1_tail_")+std::to_string(tail_ratio),"ok"};d.result.cutoff=cutoff;d.result.reference=c.reference;std::mt19937_64 rng(seed);Map current,next;current[c.observable]=1;int rz_count=0;double retained_tail_fraction_sum=0;auto start=Clock::now();
 for(auto it=c.gates.rbegin();it!=c.gates.rend();++it){
  auto advance_status=advance(current,next,*it,start,budget_s,support_cap);
  d.result.peak_terms=std::max<std::uint64_t>(d.result.peak_terms,next.size());
  if(advance_status!=AdvanceStatus::Ok){d.result.status=advance_status==AdvanceStatus::TimeCap?"time_cap":"support_cap";break;}
  current.swap(next);
  if(it->kind==GateKind::RZ&&++rz_count==interval){auto s=compress(current,cutoff,tail_ratio,rng,weighting);++d.truncation_events;d.max_heavy_terms=std::max<std::uint64_t>(d.max_heavy_terms,s.heavy);d.max_sampled_tail_terms=std::max<std::uint64_t>(d.max_sampled_tail_terms,s.sampled);d.peak_post_terms=std::max<std::uint64_t>(d.peak_post_terms,current.size());retained_tail_fraction_sum+=s.sampled_tail_l1_fraction;rz_count=0;}
  if(std::chrono::duration<double>(Clock::now()-start).count()>budget_s){d.result.status="time_cap";break;}
 }
 if(rz_count&&d.result.status=="ok"){auto s=compress(current,cutoff,tail_ratio,rng,weighting);++d.truncation_events;d.max_heavy_terms=std::max<std::uint64_t>(d.max_heavy_terms,s.heavy);d.max_sampled_tail_terms=std::max<std::uint64_t>(d.max_sampled_tail_terms,s.sampled);d.peak_post_terms=std::max<std::uint64_t>(d.peak_post_terms,current.size());retained_tail_fraction_sum+=s.sampled_tail_l1_fraction;}
 d.average_retained_tail_l1_fraction=d.truncation_events?retained_tail_fraction_sum/d.truncation_events:0;d.result.seconds=std::chrono::duration<double>(Clock::now()-start).count();d.result.memory_mb=d.result.peak_terms*48.0/1048576.0;if(d.result.status=="ok"){d.result.estimate=finish(current);d.result.error=std::abs(d.result.estimate-d.result.reference);}else d.result.error=std::numeric_limits<double>::quiet_NaN();return d;
}
}
