"""Sampling protocol (same as Colorectal): <=65 matched abstracts -> all; else 40 random (seed 42) + targeted high-risk. Random and targeted flagged separately."""
import sys, random, re, pandas as pd
key=sys.argv[1]
d=pd.read_csv(f"{key}_pool_gene_matches.csv",dtype=str).fillna(""); d["pool_rank"]=d.pool_rank.astype(int)
# context regexes: for ambiguous symbols, a document is high-risk when the gene-specific context is ABSENT (miss) or an ordinary-word context is PRESENT (hit)
MISS={"MET":r"c-?met|met (exon|amplif|mutat|inhibitor|gene|protein|receptor|proto|tyrosine|overexpress|alteration|kinase)|\bhgf\b|mesenchymal[- ]epithelial|hepatocyte growth",
 "KIT":r"c-?kit|kit (mutation|gene|receptor|exon|inhibitor|expression|d816|protein|signal|kinase|ligand)|\bgist\b|stem cell factor|cd117|imatinib|mastocyt",
 "AR":r"androgen",
 "FH":r"fumarate|hlrcc|leiomyomatosis",
 "ATM":r"ataxia|atm (mutation|gene|kinase|variant|protein|deficien|loss|inhibitor|signal|alteration)|atm-|dna (damage|repair)",
 "RET":r"ret (proto|fusion|mutation|rearrang|gene|kinase|inhibitor|receptor|alteration|ptc|mtc)|ret/ptc|medullary|selpercatinib|pralsetinib|men ?2",
 "BCR":r"bcr[- /]?abl|breakpoint cluster|philadelphia|t\(9;22\)|ph\+",
 "POLE":r"polymerase|exonuclease|pole (mutation|variant|exo|p286r|v411l|proofread)|proofreading",
 "APC":r"adenomatous|polyposis|wnt|beta-catenin|apc (mutation|gene|loss|variant|allele)|\bfap\b|catenin",
 "SET":r"set (domain|gene)|setd2",
 "CD274":r"cd274",
 "RB1":r"retinoblastoma|rb1",
 "IGH":r"immunoglobulin|igh[ /]|translocation|t\(",
 "DIS3":r"dis3","SPOP":r"spop","MYC":r"c-?myc|myc (amplif|gene|oncogene|translocation|rearrang|expression)|myc",
 "HFE":r"hemochromatosis|hfe (gene|mutation|c282y|h63d)|c282y|h63d|iron"}
HIT={"TERT":r"tert[- ]?(butyl|butoxy|amyl|butanol|butyloxy|octyl|pentyl)|tertiary butyl","APC":r"antigen[- ]presenting|activated protein c|argon plasma|\bapcs\b|dendritic cell|apc[- ]?cells?",
     "EGFR":r"(epidermal growth factor receptor|egfr)[ -]?(2\b|ii\b)",
     "BRCA2":r"partner and localizer of brca2","KIT":r"(assay|test|pcr|commercial|qiagen|elisa|extraction|diagnostic|screening|reagent|detection|manufacturer|using the|cobas|therascreen|idylla)[ -]kit|kit (was|were) (used|purchased|from)|kits?\b.{0,20}(manufacturer|purchased)"}
rng=random.Random(42); out=[]; counts=[]
for g,x in d.groupby("gene",sort=False):
    x=x.copy(); sym=g.lower()
    x["alias_only"]=x.matched_terms.apply(lambda t: sym not in t.split("|"))
    text=(x.title+" "+x.abstract).str.replace("_"," ")
    flag=pd.Series(False,index=x.index); why=pd.Series("",index=x.index)
    if g in MISS: m=~text.str.contains(MISS[g],regex=True); flag|=m; why[m]+="no-gene-context;"
    if g in HIT: m=text.str.contains(HIT[g],regex=True); flag|=m; why[m]+="ordinary-or-other-gene-context;"
    n_alias=int(x.alias_only.sum()); x["tgt_alias"]=x.alias_only; x["tgt_risk"]=flag; x["tgt_reason"]=why+x.alias_only.map({True:"alias-only;",False:""})
    # cap each targeted stratum at 40 (seeded) so large strata do not swamp the sample; the total stratum size is recorded
    sel_alias=set(x.index[x.alias_only]); sel_risk=set(x.index[flag])
    if len(sel_alias)>40: sel_alias=set(rng.sample(sorted(sel_alias),40))
    if len(sel_risk)>40: sel_risk=set(rng.sample(sorted(sel_risk),40))
    tgt=sel_alias|sel_risk
    if len(x)<=65: ran=set(); method="all"; pick=set(x.index)
    else: ran=set(rng.sample(sorted(x.index),40)); pick=ran|tgt
    s=x.loc[sorted(pick)].copy(); s["in_random"]=s.index.isin(ran); s["in_targeted"]=s.index.isin(tgt)
    s["method"]=s.apply(lambda r:"all" if len(x)<=65 else "+".join(m for m,f in (("random",r.in_random),("targeted",r.in_targeted)) if f),axis=1)
    s["n_matched"]=len(x); s["n_alias_only"]=n_alias; s["n_risk_stratum"]=len(flag[flag]); out.append(s)
    counts.append(dict(gene=g,matched=len(x),read=len(s),random=int(s.in_random.sum()),targeted=int(s.in_targeted.sum()),alias_only_total=n_alias,risk_total=int(flag.sum())))
pd.concat(out).to_csv(f"{key}_reading_set.csv",index=False,encoding="utf-8")
print(pd.DataFrame(counts).to_string(index=False)); print("total read:",sum(c["read"] for c in counts))
