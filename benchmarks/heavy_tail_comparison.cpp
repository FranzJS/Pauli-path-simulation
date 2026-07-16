#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include "pauli_bench/v1.hpp"

namespace pbv1 {
using Map = std::unordered_map<Pauli,double,PauliHash>;
using Clock = std::chrono::steady_clock;

struct RunResult {
 std::string family, case_name, method, status;
 double tail_ratio{0}, seconds{0}, memory_mb{0}, estimate{0}, reference{0}, error{0}, cutoff{0};
 std::uint64_t peak_pre_terms{0}, peak_post_terms{0}, max_heavy_k{0}, max_tail_samples{0}, truncation_events{0};
};

static double elapsed(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
static int label(const Pauli&p,int q){return (((p.x>>q)&1)?1:0)|(((p.z>>q)&1)?2:0);}
static void set_label(Pauli&p,int q,int l){auto m=1ULL<<q;p.x=(p.x&~m)|((l&1)?m:0);p.z=(p.z&~m)|((l&2)?m:0);}
static void apply_h(Pauli&p,double&w,int q){int l=label(p,q);if(l==1)set_label(p,q,2);else if(l==2)set_label(p,q,1);else if(l==3)w=-w;}
static void apply_s(Pauli&p,double&w,int q){int l=label(p,q);if(l==1){set_label(p,q,3);w=-w;}else if(l==3)set_label(p,q,1);}
struct CnotEntry{int a=0,b=0,s=1;};
static const CnotEntry& cnot_entry(int a,int b){
 static CnotEntry table[4][4]; static bool init=false;
 if(!init){using C=std::complex<double>;C I[2][2]={{1,0},{0,1}},X[2][2]={{0,1},{1,0}},Z[2][2]={{1,0},{0,-1}},Y[2][2]={{0,C(0,-1)},{C(0,1),0}};C(*P[4])[2]={I,X,Z,Y};C U[4][4]={{1,0,0,0},{0,1,0,0},{0,0,0,1},{0,0,1,0}};for(int a0=0;a0<4;++a0)for(int b0=0;b0<4;++b0){C A[4][4]{};for(int i=0;i<2;++i)for(int j=0;j<2;++j)for(int k=0;k<2;++k)for(int l=0;l<2;++l)A[2*i+k][2*j+l]=P[a0][i][j]*P[b0][k][l];C B[4][4]{};for(int i=0;i<4;++i)for(int j=0;j<4;++j)for(int k=0;k<4;++k)for(int m=0;m<4;++m)B[i][j]+=std::conj(U[k][i])*A[k][m]*U[m][j];for(int ao=0;ao<4;++ao)for(int bo=0;bo<4;++bo){C tr=0.;for(int i=0;i<2;++i)for(int j=0;j<2;++j)for(int k=0;k<2;++k)for(int l=0;l<2;++l)tr+=std::conj(P[ao][i][j]*P[bo][k][l])*B[2*i+k][2*j+l];tr/=4.;if(std::abs(std::abs(tr.real())-1)<1e-9)table[a0][b0]={ao,bo,tr.real()>0?1:-1};}}init=true;}
 return table[a][b];
}
static void apply_cnot(Pauli&p,double&w,int c,int t){auto e=cnot_entry(label(p,c),label(p,t));set_label(p,c,e.a);set_label(p,t,e.b);w*=e.s;}

static void advance(const Map&cur,Map&next,const Gate&g){
 next.clear();
 next.reserve((g.kind==GateKind::RZ?2:1)*cur.size()+1);
 for(const auto&[p0,w0]:cur){
  if(g.kind==GateKind::RZ){
   int l=label(p0,g.q0);
   if(l==1||l==3){double co=std::cos(g.theta),si=std::sin(g.theta);Pauli p2=p0;next[p0]+=w0*co;if(l==1){set_label(p2,g.q0,3);next[p2]+=w0*(-si);}else{set_label(p2,g.q0,1);next[p2]+=w0*si;}}
   else next[p0]+=w0;
  } else if(g.kind==GateKind::DEPOL){double w=w0;if(label(p0,g.q0)!=0)w*=1.0-g.theta;next[p0]+=w;}
  else {Pauli p=p0;double w=w0;if(g.kind==GateKind::H)apply_h(p,w,g.q0);else if(g.kind==GateKind::S)apply_s(p,w,g.q0);else apply_cnot(p,w,g.q0,g.q1);next[p]+=w;}
 }
}

struct SortedEntry { Pauli p; double w; double a; };
struct Split { std::vector<SortedEntry> entries; std::size_t cut{0}; };
static Split l1_split(const Map&m,double frac){
 Split s; s.entries.reserve(m.size()); long double l1=0;
 for(const auto&[p,w]:m){double a=std::abs(w); if(a>0){s.entries.push_back({p,w,a}); l1+=a;}}
 std::sort(s.entries.begin(),s.entries.end(),[](const auto&a,const auto&b){return a.a<b.a;});
 long double budget=frac*l1, used=0;
 while(s.cut<s.entries.size() && used+s.entries[s.cut].a<=budget){used+=s.entries[s.cut].a;++s.cut;}
 return s;
}

static void deterministic_truncate(Map&m,double frac,std::uint64_t&heavy_k){
 if(frac<=0||m.empty()){heavy_k=m.size();return;}
 auto s=l1_split(m,frac); heavy_k=s.entries.size()-s.cut;
 if(!s.cut)return;
 Map out; out.reserve((s.entries.size()-s.cut)*2+1);
 for(std::size_t i=s.cut;i<s.entries.size();++i)out.emplace(s.entries[i].p,s.entries[i].w);
 m.swap(out);
}

// Fixed-size PPS systematic sampling. The inclusion probabilities are
// pi_i=min(1,a_i/tau), sum pi_i=M. Randomizing order removes dependence on
// the deterministic magnitude ordering. Selected coordinates retain their
// original coefficients; this is a biased randomized truncation heuristic.
static std::size_t sample_tail_pps(const std::vector<SortedEntry>&v,std::size_t cut,std::size_t M,std::mt19937_64&rng,Map&out){
 if(cut==0||M==0)return 0;
 M=std::min(M,cut);
 long double remaining_sum=0; for(std::size_t i=0;i<cut;++i)remaining_sum+=v[i].a;
 std::size_t saturated=0;
 long double tau=0;
 while(true){
  std::size_t unsat=cut-saturated;
  std::size_t slots=M-saturated;
  if(slots==0){tau=std::numeric_limits<long double>::infinity();break;}
  tau=remaining_sum/static_cast<long double>(slots);
  if(unsat==0 || static_cast<long double>(v[unsat-1].a)<=tau*(1.0L+1e-18L))break;
  remaining_sum-=v[unsat-1].a; ++saturated;
 }
 std::vector<std::size_t> order(cut); for(std::size_t i=0;i<cut;++i)order[i]=i; std::shuffle(order.begin(),order.end(),rng);
 std::vector<long double> pi(cut);
 long double sum_pi=0;
 for(std::size_t i=0;i<cut;++i){pi[i]=std::min<long double>(1.0L,static_cast<long double>(v[i].a)/tau);sum_pi+=pi[i];}
 long double drift=static_cast<long double>(M)-sum_pi;
 for(std::size_t k=0;k<cut && std::abs(drift)>1e-14L;++k){auto i=order[k];long double room=drift>0?1.0L-pi[i]:pi[i];long double d=std::copysign(std::min(std::abs(drift),room),drift);pi[i]+=d;drift-=d;}
 std::uniform_real_distribution<double> u(0.0,1.0); long double point=u(rng), cumulative=0; std::size_t selected=0;
 for(auto idx:order){long double next=cumulative+pi[idx];if(selected<M && point<next){out.emplace(v[idx].p,v[idx].w);++selected;point+=1.0L;}cumulative=next;}
 if(selected!=M){
  std::vector<std::size_t> candidates(order);std::sort(candidates.begin(),candidates.end(),[&](auto a,auto b){return pi[a]>pi[b];});
  for(auto idx:candidates){if(selected>=M)break;if(out.find(v[idx].p)==out.end()){out.emplace(v[idx].p,v[idx].w);++selected;}}
 }
 return selected;
}

static void hybrid_truncate_fixed_k(Map&m,std::size_t target_k,double ratio,std::mt19937_64&rng,std::uint64_t&heavy_k,std::uint64_t&tail_samples){
 if(m.empty()){heavy_k=tail_samples=0;return;}
 std::vector<SortedEntry> entries;entries.reserve(m.size());
 for(const auto&[p,w]:m){double a=std::abs(w);if(a>0)entries.push_back({p,w,a});}
 std::sort(entries.begin(),entries.end(),[](const auto&a,const auto&b){return a.a<b.a;});
 const std::size_t K=std::min(target_k,entries.size());const std::size_t cut=entries.size()-K;heavy_k=K;
 if(cut==0||K==0){tail_samples=0;return;}
 const std::size_t M=std::min<std::size_t>(cut,std::max<std::size_t>(1,static_cast<std::size_t>(std::ceil(ratio*K))));
 Map out;out.reserve((K+M)*2+1);for(std::size_t i=cut;i<entries.size();++i)out.emplace(entries[i].p,entries[i].w);
 tail_samples=sample_tail_pps(entries,cut,M,rng,out);m.swap(out);
}

static double finish(const Map&m){long double v=0;for(const auto&[p,w]:m)if(p.x==0)v+=w;return static_cast<double>(v);}

static RunResult run(const Circuit&c,double budget,double ratio,std::uint64_t seed,std::vector<std::size_t>*schedule_out=nullptr,const std::vector<std::size_t>*schedule_in=nullptr){
 RunResult r{c.family,c.name,ratio==0?"bfs_l1":"bfs_l1_heavy_plus_random_tail","ok",ratio};r.reference=c.reference;r.cutoff=c.cutoff;
 auto st=Clock::now();std::mt19937_64 rng(seed);Map cur,next;cur[c.observable]=1.;r.peak_pre_terms=r.peak_post_terms=1;int count=0;std::size_t event=0;
 auto compress=[&](){std::uint64_t K=0,M=0;if(ratio==0){deterministic_truncate(cur,c.cutoff,K);if(schedule_out)schedule_out->push_back(K);}else{if(!schedule_in||event>=schedule_in->size()){r.status="schedule_error";return;}hybrid_truncate_fixed_k(cur,(*schedule_in)[event],ratio,rng,K,M);}r.max_heavy_k=std::max(r.max_heavy_k,K);r.max_tail_samples=std::max(r.max_tail_samples,M);r.peak_post_terms=std::max<std::uint64_t>(r.peak_post_terms,cur.size());++r.truncation_events;++event;};
 for(auto it=c.gates.rbegin();it!=c.gates.rend();++it){advance(cur,next,*it);cur.swap(next);r.peak_pre_terms=std::max<std::uint64_t>(r.peak_pre_terms,cur.size());if(it->kind==GateKind::RZ&&++count==4){compress();count=0;if(r.status!="ok")break;}if(elapsed(st)>budget){r.status="time_cap";break;}}
 if(count&&r.status=="ok")compress();
 if(r.status=="ok")r.estimate=finish(cur);
 r.seconds=elapsed(st);r.memory_mb=r.peak_pre_terms*48.0/1048576.0;r.error=r.status=="ok"?std::abs(r.estimate-r.reference):std::numeric_limits<double>::quiet_NaN();return r;
}

}

int main(int argc,char**argv){
 using namespace pbv1; int only=argc>1?std::stoi(argv[1]):-1; int passes=argc>2?std::stoi(argv[2]):10; double budget=argc>3?std::stod(argv[3]):180.0;
 std::ofstream one("results/hybrid_single_pass.csv");
 std::ofstream prog("results/hybrid_progressive_10_pass.csv");
 one<<"family,case,method,tail_ratio,pass,seed,status,runtime_s,frontier_memory_mb,peak_pre_terms,peak_post_terms,max_heavy_k,max_tail_samples,truncation_events,estimate,reference,abs_error,cutoff\n";
 prog<<"family,case,tail_ratio,passes,mean_estimate,abs_error,empirical_sd,empirical_stderr,mean_runtime_s,total_runtime_s,max_frontier_memory_mb,reference\n";
 auto circuits=benchmark_circuits(); const std::vector<double> ratios{0.05,0.10,0.20};
 for(int ci=0;ci<(int)circuits.size();++ci){if(only>=0&&ci!=only)continue;auto&c=circuits[ci];
  std::vector<std::size_t> schedule;auto base=run(c,budget,0,20260716,&schedule,nullptr); one<<base.family<<','<<base.case_name<<','<<base.method<<",0,1,20260716,"<<base.status<<','<<std::setprecision(17)<<base.seconds<<','<<base.memory_mb<<','<<base.peak_pre_terms<<','<<base.peak_post_terms<<','<<base.max_heavy_k<<','<<base.max_tail_samples<<','<<base.truncation_events<<','<<base.estimate<<','<<base.reference<<','<<base.error<<','<<base.cutoff<<'\n';one.flush();
  std::cout<<"BASE "<<c.name<<" t="<<base.seconds<<" mem="<<base.memory_mb<<" err="<<base.error<<"\n"<<std::flush;
  for(double ratio:ratios){long double sum=0,sum2=0,total_time=0;double maxmem=0;for(int p=1;p<=passes;++p){std::uint64_t seed=202607160000ULL+1009ULL*ci+1000003ULL*static_cast<std::uint64_t>(std::llround(ratio*100))+p;auto x=run(c,budget,ratio,seed,nullptr,&schedule);if(p==1){one<<x.family<<','<<x.case_name<<','<<x.method<<','<<ratio<<','<<p<<','<<seed<<','<<x.status<<','<<std::setprecision(17)<<x.seconds<<','<<x.memory_mb<<','<<x.peak_pre_terms<<','<<x.peak_post_terms<<','<<x.max_heavy_k<<','<<x.max_tail_samples<<','<<x.truncation_events<<','<<x.estimate<<','<<x.reference<<','<<x.error<<','<<x.cutoff<<'\n';one.flush();}
    long double oldmean=p>1?sum/(p-1):0;sum+=x.estimate;long double newmean=sum/p;if(p>1)sum2+=(x.estimate-oldmean)*(x.estimate-newmean);total_time+=x.seconds;maxmem=std::max(maxmem,x.memory_mb);double sd=p>1?std::sqrt(static_cast<double>(sum2/(p-1))):std::numeric_limits<double>::quiet_NaN();double se=p>1?sd/std::sqrt((double)p):std::numeric_limits<double>::quiet_NaN();prog<<c.family<<','<<c.name<<','<<ratio<<','<<p<<','<<std::setprecision(17)<<static_cast<double>(newmean)<<','<<std::abs(static_cast<double>(newmean)-c.reference)<<','<<sd<<','<<se<<','<<static_cast<double>(total_time/p)<<','<<static_cast<double>(total_time)<<','<<maxmem<<','<<c.reference<<'\n';prog.flush();std::cout<<"PASS "<<c.name<<" r="<<ratio<<" p="<<p<<" t="<<x.seconds<<" est="<<x.estimate<<" avg_err="<<std::abs(static_cast<double>(newmean)-c.reference)<<"\n"<<std::flush;}
  }
 }
}
