"""All-19 synthesis of the top-10 gene-mention QC. Aggregates <key>_qc_final.csv (reader identity + evidence). Pool-wide numbers come only from
complete alias strata that were read in full, or from the stated automated scans; nothing here is extrapolated from the random samples."""
import pandas as pd, re, glob
KEYS=[("breast","Breast Cancer"),("lung","Lung Cancer"),("colorectal","Colorectal Cancer"),("prostate","Prostate Cancer"),("melanoma","Melanoma"),("bladder","Bladder Cancer"),
 ("kidney","Kidney Cancer"),("pancreas","Pancreatic Cancer"),("liver","Liver Cancer"),("stomach","Stomach Cancer"),("esophageal","Esophageal Cancer"),("ovarian","Ovarian Cancer"),
 ("endometrial","Endometrial Cancer"),("cervical","Cervical Cancer"),("thyroid","Thyroid Cancer"),("brain","Brain Cancer"),("leukemia","Leukemia"),("lymphoma","Lymphoma"),("myeloma","Myeloma")]
frames=[];pools=[]
for k,l in KEYS:
    f=pd.read_csv(f"{k}_qc_final.csv",dtype=str).fillna(""); f["key"]=k; frames.append(f)
    p=pd.read_csv(f"{k}_pool_gene_matches.csv",dtype=str).fillna(""); p["key"]=k; p["cancer"]=l; pools.append(p[["key","cancer","gene","pmid"]])
F=pd.concat(frames,ignore_index=True); P=pd.concat(pools,ignore_index=True)
F["protocol"]=F.key.map(lambda k:"v1 (top15 + random25 + alias)" if k in("breast","lung") else "v2 (<=65 all; else 40 random + targeted)")
F["stratum"]=F.sampling_method.map(lambda m:"all-matched" if m=="all" else ("random+targeted" if "random" in m and "targeted" in m else ("random" if "random" in m else ("targeted-only" if "targeted" in m else ("top15(v1)" if "top15" in m else "alias(v1)")))))
credits=P.groupby(["cancer","gene"]).size().rename("pool_credits").reset_index()
# ---- per-cancer
def cs(g):
    return pd.Series(dict(top10_genes=g.gene.nunique(),abstracts_read=len(g),random_or_all=int(g.stratum.isin(["random","random+targeted","all-matched"]).sum()),targeted_only=int((g.stratum=="targeted-only").sum()),
        incorrect=int((g.identity=="X").sum()),uncertain=int((g.identity=="U").sum())))
C=F.groupby("cancer",sort=False).apply(cs).reset_index(); C=C.merge(P.groupby("cancer",sort=False).size().rename("pool_credits_top10").reset_index(),on="cancer"); C.to_csv("ALL19_cancer_summary.csv",index=False)
# ---- mechanisms
def mech(r):
    t=r.rationale.lower()
    for key,lab in (("'pi3k'","PI3K -> PIK3CA (pathway/family alias)"),("hnpcc","HNPCC syndrome name -> MLH1/MSH2"),("inside her2","EGFR credited through HER2's name"),("receptor' inside her2","EGFR credited through HER2's name"),
        ("hif-1","HIF-1 -> SETD2 (hypoxia-inducible factor)"),("wagr","WAGR syndrome name -> WT1"),("family/drug-class","PARP (family/drug-class term) -> PARP1"),("kit'","'kit' = assay kit -> KIT"),("tert-butyl","'tert' = tert-butyl chemistry -> TERT"),
        ("apc/c","APC/C (anaphase-promoting complex) -> APC"),("b-cell","BCR = B-cell receptor -> BCR"),("methionine","MET = methionine -> MET"),("metastasis","MET = metastasis -> MET"),("inside palb2","symbol inside another gene's full name (PALB2/BRCA)"),
        ("inside another gene","symbol inside another gene's full name (PALB2/BRCA)"),("men1","MEN1 used as syndrome label"),("armadillo","'armadillo' -> CTNNB1")):
        if key in t: return lab
    return "other: "+r.rationale[:60]
X=F[F.identity.isin(["X","U"])].copy(); X["mechanism"]=X.apply(mech,axis=1)
X["cancer_gene"]=X.cancer+" / "+X.gene
M=X.groupby("mechanism").apply(lambda g:pd.Series(dict(abstracts_wrong_or_uncertain=len(g),wrong=int((g.identity=="X").sum()),uncertain=int((g.identity=="U").sum()),
    cancers=len(g.cancer.unique()),detail="; ".join(f"{c}/{ge} {n}" for (c,ge),n in g.groupby(["cancer","gene"]).size().items())))).reset_index().sort_values("abstracts_wrong_or_uncertain",ascending=False)
# denominators: pool credits of affected (cancer,gene)
def denom(row):
    return ""
M.to_csv("ALL19_error_mechanisms.csv",index=False)
X.drop(columns=["cancer_gene"]).to_csv("ALL19_suspected_false_matches.csv",index=False)
# ---- affected gene credit rates
A=X.groupby(["cancer","gene","mechanism"]).size().rename("n").reset_index().merge(credits,on=["cancer","gene"]); A["pct_of_gene_credits"]=(100*A.n/A.pool_credits).round(1)
A.sort_values("pct_of_gene_credits",ascending=False).to_csv("ALL19_affected_gene_rates.csv",index=False)
# ---- evidence distribution (random / all only; targeted kept separate)
EV=["M","ML","G","P","W","I"]
def evd(g): 
    n=len(g); return pd.Series({**{e:int((g.evidence_final==e).sum()) for e in EV},"n":n,"target_ctx":int((g.target_context=="target").sum())})
R=F[F.stratum.isin(["random","random+targeted","all-matched"])]
T=F[F.stratum=="targeted-only"]
E_r=R.groupby(["cancer","gene"],sort=False).apply(evd).reset_index(); E_r["pct_mutation_specific(M)"]=(100*E_r.M/E_r.n).round(0); E_r["pct_mut_or_listed(M+ML)"]=(100*(E_r.M+E_r.ML)/E_r.n).round(0)
E_r.merge(credits,on=["cancer","gene"]).to_csv("ALL19_evidence_random_or_all.csv",index=False)
E_t=T.groupby(["cancer","gene"],sort=False).apply(evd).reset_index(); E_t.to_csv("ALL19_evidence_targeted_only.csv",index=False)
tot=R.groupby("evidence_final").size().reindex(EV).fillna(0).astype(int)
print("ABSTRACTS READ:",len(F),"| random/all:",len(R),"| targeted-only:",len(T),"| top15(v1):",int((F.stratum=='top15(v1)').sum()))
print("identity (all read):",F.identity.value_counts().to_dict())
print("identity (random/all only):",R.identity.value_counts().to_dict())
print("identity (targeted-only):",T.identity.value_counts().to_dict())
print("evidence classes (random/all):",tot.to_dict(), "target_ctx:",int((R.target_context=='target').sum()))
print(); print(C.to_string(index=False)); print(); print(M[["mechanism","abstracts_wrong_or_uncertain","wrong","uncertain","cancers"]].to_string(index=False))
low=E_r[(E_r.n>=20)].sort_values("pct_mutation_specific(M)").head(25); print(); print(low[["cancer","gene","n","M","ML","G","P","W","I","pct_mutation_specific(M)","target_ctx"]].to_string(index=False))
