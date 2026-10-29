"""Analysis of the read-only counterfactuals. Compares each variant with the unmodified baseline of the SAME harness, verifies each baseline against the frozen v5
outputs, and scores the benchmark (IntOGen 'published' reference: any-cancer P@10 and exact concordance) by substituting the cancer's top-10 into the canonical 190.
Reports every changed position; no recommendation is made from the benchmark."""
import os
import json, sys, math
from pathlib import Path
import pandas as pd
REPO=Path(os.environ.get("ONCODIGGER_PAPER_REPO") or Path(__file__).resolve().parents[3]); sys.path.insert(0,str(REPO/"scripts/analysis"))
import intogen_statistical_support_rescoring as base
anyg,exact,_=base.load_references()["published"]
LAB={"breast":"Breast Cancer","lung":"Lung Cancer","colorectal":"Colorectal Cancer","prostate":"Prostate Cancer","bladder":"Bladder Cancer","kidney":"Kidney Cancer","liver":"Liver Cancer","stomach":"Stomach Cancer","esophageal":"Esophageal Cancer","endometrial":"Endometrial Cancer","cervical":"Cervical Cancer","thyroid":"Thyroid Cancer"}
canon=pd.read_csv(REPO/"data/canonical/precision_at_k_detail.csv")[["cancer","rank","gene"]]
frozen=pd.read_csv(REPO/"data/derived/query_planes/mutation_vocab_corrected.csv")
exb=(pd.read_csv(REPO/"data/derived/intogen_rescoring/intogen_literature_weight_audit.csv").dropna(subset=["engine_rank"]).rename(columns={"engine_rank":"rank","supporting_papers":"S","corpus_papers":"B_exact","relevance_score":"relevance_new"}))
def bench(c,top10):
    t=canon[(canon["rank"]<=10)&(canon.cancer!=LAB[c])]; add=pd.DataFrame({"cancer":LAB[c],"rank":range(1,11),"gene":top10}); t=pd.concat([t,add])
    t["gene"]=t.gene.str.upper(); a=int(t.gene.isin(anyg).sum()); e=int(sum((g,cc) in exact for g,cc in zip(t.gene,t.cancer)))
    cc=t[t.cancer==LAB[c]]; return a,e,int(cc.gene.isin(anyg).sum()),int(sum((g,LAB[c]) in exact for g in cc.gene))
def load(c,v):
    f=Path(f"cf_{c}_{v}.json"); 
    if c=="breast" and v=="baseline" and not f.exists(): f=Path("cf_baseline.json")
    return json.load(open(f)) if f.exists() else None
def E(o): N=o.get("N",1000); C=o.get("C",250061); return (o["support"]/N)/(o["background"]/C)
VARS={"no_pi3k":("PIK3CA","'PI3K' alias removed from PIK3CA","confirmed false identity"),"no_hif1":("SETD2","'HIF-1' alias removed from SETD2","confirmed false identity"),
      "no_hnpcc":("MLH1/MSH2","'HNPCC' syndrome alias removed from MLH1 and MSH2","confirmed false identity"),"no_parp":("PARP1","family alias 'PARP' removed from PARP1","broader family-alias policy test (NOT a confirmed-false correction)")}
ver=[];rows=[];moves=[]
for c in LAB:
    b=load(c,"baseline")
    if b is None: continue
    fz=frozen[frozen.cancer==LAB[c]].sort_values("rank").gene.tolist(); bg=[o["gene"] for o in b]
    ex=exb[exb.cancer==LAB[c]].sort_values("rank").set_index("gene")
    same_order=bg[:50]==fz[:50]; n_ord=sum(1 for x,y in zip(bg,fz) if x==y)
    sc=[(o["support"]==ex.loc[o["gene"],"S"], o["background"]==ex.loc[o["gene"],"B_exact"], abs(o["relevance"]-ex.loc[o["gene"],"relevance_new"])<1e-3) for o in b if o["gene"] in ex.index]
    m=bench(c,bg[:10])
    ver.append(dict(cancer=LAB[c],order_top50_identical_to_frozen=same_order,positions_matching=f"{n_ord}/50",support_equal=f"{sum(s[0] for s in sc)}/{len(sc)}",background_equal=f"{sum(s[1] for s in sc)}/{len(sc)}",relevance_equal_tol=f"{sum(s[2] for s in sc)}/{len(sc)}",
        top10_equals_canonical=bg[:10]==canon[(canon.cancer==LAB[c])&(canon['rank']<=10)].sort_values('rank').gene.tolist(),intogen_any_190=m[0],intogen_exact_190=m[1]))
    for v,(gene,desc,kind) in VARS.items():
        a=load(c,v)
        if a is None: continue
        pos_b={o["gene"]:o for o in b}; pos_a={o["gene"]:o for o in a}
        t10b=[o["gene"] for o in b[:10]]; t10a=[o["gene"] for o in a[:10]]
        ma=bench(c,t10a)
        out=[g for g in t10b if g not in t10a]; inn=[g for g in t10a if g not in t10b]
        genes=gene.split("/") 
        for g in genes:
            ob=pos_b.get(g); oa=pos_a.get(g)
            rows.append(dict(cancer=LAB[c],scenario=desc,kind=kind,gene=g,rank_before=ob["rank"] if ob else None,rank_after=oa["rank"] if oa else None,
              support_before=ob["support"] if ob else None,support_after=oa["support"] if oa else None,background_before=ob["background"] if ob else None,background_after=oa["background"] if oa else None,
              enrich_before=round(E(ob),3) if ob else None,enrich_after=round(E(oa),3) if oa else None,relevance_before=ob["relevance"] if ob else None,relevance_after=oa["relevance"] if oa else None,
              top10_replacements=("; ".join(f"{o} -> {i}" for o,i in zip(out,inn)) if out else "none"),
              intogen_any_before_after=f"{bench(c,t10b)[0]} -> {ma[0]}",intogen_exact_before_after=f"{bench(c,t10b)[1]} -> {ma[1]}",
              cancer_any_exact_before_after=f"{bench(c,t10b)[2]}/{bench(c,t10b)[3]} -> {ma[2]}/{ma[3]} (of 10)"))
        for g in sorted(set(pos_b)|set(pos_a)):
            rb=pos_b[g]["rank"] if g in pos_b else None; ra=pos_a[g]["rank"] if g in pos_a else None
            if rb!=ra: moves.append(dict(cancer=LAB[c],scenario=desc,gene=g,rank_before=rb,rank_after=ra))
pd.DataFrame(ver).to_csv("CF_baseline_verification.csv",index=False); R=pd.DataFrame(rows); R.to_csv("CF_consolidated.csv",index=False); pd.DataFrame(moves).to_csv("CF_all_position_changes.csv",index=False)
pd.set_option("display.width",250); pd.set_option("display.max_colwidth",60)
print("BASELINE VERIFICATION"); print(pd.DataFrame(ver).to_string(index=False)); print(); print("CONSOLIDATED"); print(R.drop(columns=["kind"]).to_string(index=False))
