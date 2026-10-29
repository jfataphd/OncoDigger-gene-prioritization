"""Merge reader identity verdicts with the evidence draft. usage: finalize.py <key> <Cancer label>
Identity rules were verified by reading every flagged abstract (see <key> notes); default identity = C (correct) when the symbol or an unambiguous
full name appears. X = incorrect, U = uncertain."""
import sys, re, pandas as pd
key,label=sys.argv[1],sys.argv[2]
s=pd.read_csv(f"{key}_reading_set_classified.csv",dtype=str).fillna("")
ov=pd.read_csv(f"{key}_verdict_overrides.csv",dtype=str).fillna("") if __import__("os").path.exists(f"{key}_verdict_overrides.csv") else pd.DataFrame(columns=["gene","pmid","identity","evidence","rationale"])
VAR={"MLH1":r"h?mlh[ -]?1|mutl homolog 1","MSH2":r"h?msh[ -]?2|muts homolog 2","MSH6":r"h?msh[ -]?6|muts homolog 6"}
rows=[]
for _,r in s.iterrows():
    text=(r.title+" "+r.abstract).replace("_"," "); sym=r.gene.lower(); terms=r.matched_terms.split("|")
    alias_only=sym not in terms
    ident,why="C",f"symbol '{r.gene}' present" if not alias_only else f"credited via alias {'|'.join(terms)}"
    if alias_only:
        if r.gene in VAR and "hnpcc" in terms:
            if re.search(r"(?<![a-z0-9])(%s)(?![a-z0-9])"%VAR[r.gene],text): ident,why="C","credited via 'HNPCC' but a variant spelling of the gene (e.g. hMLH1/hMSH2) is in the text"
            else: ident,why="X","credited only via the syndrome name 'HNPCC'; no spelling of the gene in the abstract"
        elif r.gene=="PIK3CA" and terms==["pi3k"]: ident,why="X","credited only via the alias 'PI3K' (pathway/family term); PIK3CA not named"
        elif r.gene=="EGFR" and re.search(r"receptor[ -]?2\b",text): ident,why="X","credited via 'epidermal growth factor receptor' inside HER2's name"
        elif r.gene=="SETD2" and terms==["hif 1"]: ident,why="X","credited only via the alias 'HIF-1', which in these abstracts is hypoxia-inducible factor 1"
        elif r.gene=="WT1" and terms==["wagr"]: ident,why="X","credited only via the syndrome name 'WAGR'; the abstract is about PAX6/aniridia, WT1 not named"
        elif r.gene=="PARP1" and terms==["parp"]:
            if re.search(r"parp[ -]?1",text): ident,why="C","'PARP-1' written in the text"
            else: ident,why="U","credited only via the family/drug-class term 'PARP' (e.g. 'PARP inhibitor'); PARP1 not specified"
        elif r.gene=="ERBB2": ident,why="C","HER2/HER-2 = ERBB2 product (alias); see evidence class for subtype/status usage"
        elif r.gene=="CD274": ident,why="C","PD-L1 = CD274 product (alias)"
    rows.append(dict(cancer=label,gene=r.gene,pmid=r.pmid,sampling_method=r.method,triggering_alias="|".join(terms),identity=ident,evidence_draft=r.evidence_draft,target_context=r.target_context,rationale=why,pool_rank=r.pool_rank,title=r.title[:140]))
f=pd.DataFrame(rows)
for _,o in ov.iterrows():
    m=(f.gene==o.gene)&(f.pmid==o.pmid)
    if o.identity: f.loc[m,"identity"]=o.identity
    if o.evidence: f.loc[m,"evidence_draft"]=o.evidence
    if o.rationale: f.loc[m,"rationale"]=o.rationale
f["evidence_final"]=f.evidence_draft
f.loc[(f.gene=="PARP1")&(f.identity=="U"),"evidence_final"]="W"
f.to_csv(f"{key}_qc_final.csv",index=False,encoding="utf-8")
def stratum(df): return pd.Series(dict(n=len(df),correct=(df.identity=="C").sum(),incorrect=(df.identity=="X").sum(),uncertain=(df.identity=="U").sum()))
f["is_random"]=f.sampling_method.str.contains("random|all"); f["is_targeted"]=f.sampling_method.str.contains("targeted")
out=[]
for g,x in f.groupby("gene",sort=False):
    rr=x[x.is_random]; tt=x[x.is_targeted&~x.is_random]
    out.append(dict(gene=g,read=len(x),random_n=len(rr),random_incorrect=int((rr.identity=="X").sum()),targeted_only_n=len(tt),targeted_incorrect=int((tt.identity=="X").sum()),
        uncertain=int((x.identity=="U").sum()),**{f"ev_{c}":int((x.evidence_final==c).sum()) for c in ("M","ML","G","P","W","I")},target_ctx=int((x.target_context=="target").sum())))
S=pd.DataFrame(out); S.to_csv(f"{key}_qc_gene_summary.csv",index=False); print(S.to_string(index=False))
f[f.identity!="C"].drop(columns=["is_random","is_targeted"]).to_csv(f"{key}_suspected_false_matches.csv",index=False)
print("suspected false matches:",int((f.identity!="C").sum()))
