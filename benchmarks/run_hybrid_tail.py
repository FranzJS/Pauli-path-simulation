#!/usr/bin/env python3
"""Run process-isolated BFS versus L1-heavy-hitters + randomized-tail benchmarks."""
from __future__ import annotations
import argparse,csv,math,statistics,subprocess
from collections import defaultdict
from pathlib import Path
FIELDS=["family","case","method","tail_ratio","seed","status","runtime_s","memory_mb","peak_pre_terms","peak_post_terms","max_heavy","max_tail","tail_l1_fraction","estimate","reference","error","events"]
RATIOS=(0.05,0.10,0.20)
RATIO_OFFSET={0.05:5045,0.10:10090,0.20:20180}
DEFAULT_CAP={(0,0.05):1500000,(0,0.10):1500000,(0,0.20):1500000,(1,0.05):800000,(1,0.10):800000,(1,0.20):700000,(2,0.05):1500000,(2,0.10):1500000,(2,0.20):1500000}
NUMERIC_FLOAT={"tail_ratio","runtime_s","memory_mb","tail_l1_fraction","estimate","reference","error"}
NUMERIC_INT={"seed","peak_pre_terms","peak_post_terms","max_heavy","max_tail","events"}
def parse_row(text:str)->dict:
 row=next(csv.DictReader([",".join(FIELDS),text.strip()]))
 for k in NUMERIC_FLOAT: row[k]=float(row[k]) if row[k] not in ("","nan") else math.nan
 for k in NUMERIC_INT: row[k]=int(row[k]) if row[k] else 0
 return row
def run(exe:Path,args:list[str],timeout:float)->dict:
 try:p=subprocess.run([str(exe),*args],text=True,capture_output=True,timeout=timeout,check=False)
 except subprocess.TimeoutExpired:return {**{k:"" for k in FIELDS},"method":"hybrid_l1_tail","status":"external_timeout"}
 if not p.stdout.strip():raise RuntimeError(p.stderr.strip() or "benchmark emitted no row")
 return parse_row(p.stdout.strip().splitlines()[-1])
def write_csv(path:Path,rows:list[dict],fields:list[str])->None:
 with path.open("w",newline="") as f:w=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore");w.writeheader();w.writerows(rows)
def main()->None:
 ap=argparse.ArgumentParser();ap.add_argument("--exe",type=Path,default=Path("build/pauli_hybrid_tail_one"));ap.add_argument("--output-dir",type=Path,default=Path("results"));ap.add_argument("--passes",type=int,default=10);ap.add_argument("--budget",type=float,default=120.0);ap.add_argument("--seed",type=int,default=20260716);ap.add_argument("--weighting",choices=("original","ht"),default="original");args=ap.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
 baselines=[];raw=[]
 for ci in range(3):baselines.append(run(args.exe,["bfs",str(ci),str(args.budget)],args.budget+60))
 for ci in range(3):
  for ratio in RATIOS:
   for p in range(args.passes):
    seed=args.seed+1000003*p+RATIO_OFFSET[ratio]+7919*ci;cap=DEFAULT_CAP[(ci,ratio)]
    row=run(args.exe,["hybrid",str(ci),str(ratio),str(args.budget),str(seed),str(cap),args.weighting],args.budget+60);raw.append(row)
    if row.get("status")!="ok":break
 write_csv(args.output_dir/"hybrid_l1_tail_raw.csv",raw,FIELDS)
 baseline_by_case={r["case"]:r for r in baselines};single=[]
 for b in baselines:single.append({**b,"pass":1,"runtime_ratio_vs_bfs":1.0,"memory_ratio_vs_bfs":1.0})
 grouped=defaultdict(list)
 for r in raw:grouped[(r.get("case"),r.get("tail_ratio"))].append(r)
 for group in grouped.values():
  r=group[0];b=baseline_by_case[r["case"]];single.append({**r,"pass":1,"runtime_ratio_vs_bfs":r["runtime_s"]/b["runtime_s"],"memory_ratio_vs_bfs":r["memory_mb"]/b["memory_mb"]})
 single_fields=FIELDS+["pass","runtime_ratio_vs_bfs","memory_ratio_vs_bfs"];write_csv(args.output_dir/"hybrid_l1_tail_single_pass.csv",single,single_fields)
 cumulative=[]
 for (case,ratio),group in grouped.items():
  estimates=[];runtime=0.0;peak=0.0
  for i,r in enumerate(group,1):
   runtime+=float(r.get("runtime_s") or 0);peak=max(peak,float(r.get("memory_mb") or 0))
   if r.get("status")!="ok":cumulative.append({"family":r.get("family"),"case":case,"method":r.get("method"),"tail_ratio":ratio,"passes":i,"status":r.get("status"),"cumulative_runtime_s":runtime,"mean_runtime_s":"","peak_memory_mb":peak,"cumulative_mean":"","reference":r.get("reference"),"abs_error_of_cumulative_mean":"","sample_stddev":"","standard_error":""});break
   estimates.append(r["estimate"]);mean=statistics.fmean(estimates);sd=statistics.stdev(estimates) if i>1 else math.nan
   cumulative.append({"family":r["family"],"case":case,"method":r["method"],"tail_ratio":ratio,"passes":i,"status":"ok","cumulative_runtime_s":runtime,"mean_runtime_s":runtime/i,"peak_memory_mb":peak,"cumulative_mean":mean,"reference":r["reference"],"abs_error_of_cumulative_mean":abs(mean-r["reference"]),"sample_stddev":sd,"standard_error":sd/math.sqrt(i) if i>1 else math.nan})
 cum_fields=["family","case","method","tail_ratio","passes","status","cumulative_runtime_s","mean_runtime_s","peak_memory_mb","cumulative_mean","reference","abs_error_of_cumulative_mean","sample_stddev","standard_error"];write_csv(args.output_dir/"hybrid_l1_tail_10_pass.csv",cumulative,cum_fields)
if __name__=="__main__":main()
