"""Rule-based DRAFT evidence class for each sampled abstract (reviewed/overridden by the reader). usage: classify.py <key> [reading_set.csv]
Classes (priority): M mutation evidence specific to the matched gene; ML mutation context but gene only in a multi-gene list; G other gene-specific biology
(methylation, amplification, copy number, mRNA expression, functional/inhibitor); P protein expression/status (IHC, biomarker status); W pathway/family mention
(no symbol, or 'pathway' wording); I incidental. Also flag target_context (named as the target of an antibody/inhibitor, e.g. anti-EGFR therapy)."""
import os
from pathlib import Path
import sys, re, pandas as pd
key=sys.argv[1]; f=sys.argv[2] if len(sys.argv)>2 else f"{key}_reading_set.csv"
s=pd.read_csv(f,dtype=str).fillna("")
SYMS=set(x.strip().lower() for x in pd.read_csv(str(Path(os.environ.get("ONCODIGGER_PAPER_REPO") or Path(__file__).resolve().parents[3]) / "reference_engine_code" / "Human Genes.csv"),encoding="utf-8-sig",usecols=["symbol"]).symbol if "-" not in str(x) and 3<=len(str(x))<=8)
MUT=re.compile(r"mutat|variant|alteration|deletion|insertion|substitution|frameshift|nonsense|missense|splice|truncat|hotspot|germline|somatic|loss[- ]of[- ]function|gain[- ]of[- ]function|wild[- ]type|exon \d+|\bv600|codon|point mutation|\bsnv\b|indel|pathogenic|carrier",re.I)
G=re.compile(r"methylat|promoter|mrna|transcript|amplif|copy[- ]number|fusion|rearrang|translocation|knock(down|out)|crispr|sirna|shrna|inhibitor|polymorphism|\bsnp\b|haplotype|silenc|epigen|overexpress|gene expression|signature|resistance",re.I)
P=re.compile(r"\bihc\b|immunohisto|protein|positiv|negativ|status|staining|expression|deficien|loss of|level",re.I)
W=re.compile(r"pathway|signal+ing|cascade|family|axis|network",re.I)
T=re.compile(r"anti-?(epidermal growth factor receptor|egfr)|egfr[- ](targeted|inhibitor|antibod|tki)|cetuximab|panitumumab|trastuzumab|her2[- ]targeted|anti-?her2|\btki\b|inhibitors? (of|targeting)|targeting (the )?",re.I)
def terms(r): return [t.replace("_"," ") for t in r.matched_terms.split("|") if t]
rows=[]
for _,r in s.iterrows():
    text=(r.title+" || "+r.abstract).replace("_"," "); sym=r.gene.lower(); tl=terms(r)
    spans=[m for t in tl for m in re.finditer(r"(?<![a-z0-9])"+re.escape(t).replace(r"\ ","[ -]")+r"(?![a-z0-9])",text)]
    if not spans: rows.append(("I","","no-match-in-text")); continue
    flags=set(); tgt=False; listed_all=True
    for m in spans:
        w=text[max(0,m.start()-90):m.end()+90]; tight=text[max(0,m.start()-35):m.end()+35]; toks=re.findall(r"[a-z0-9]+",w)
        ngene=sum(1 for t in toks if t in SYMS and t!=m.group(0)); inlist=ngene>=3
        t_hit=bool(T.search(w)); tgt|=t_hit
        if MUT.search(tight): flags.add("ML" if inlist else "M")
        elif t_hit: flags.add("I")
        elif G.search(tight): flags.add("G")
        elif P.search(tight): flags.add("P")
        elif W.search(tight): flags.add("W")
        else: flags.add("I")
    for c in ("M","G","P","W","ML","I"):
        if c in flags: cls=c; break
    aliasonly=sym not in [t for t in r.matched_terms.split("|")]
    if aliasonly and cls in("M","ML","G","P") and not any(re.search(r"(?<![a-z0-9])"+sym+r"(?![a-z0-9])",text) for _ in [0]) and any(t in("pi3k","hnpcc","pd l1","her2","her 2") for t in tl): cls_note="alias-context"
    else: cls_note=""
    rows.append((cls,"target" if tgt else "",cls_note))
s["evidence_draft"]=[x[0] for x in rows]; s["target_context"]=[x[1] for x in rows]; s["draft_note"]=[x[2] for x in rows]
s.to_csv(f"{key}_reading_set_classified.csv",index=False,encoding="utf-8")
print(s.groupby("gene",sort=False).evidence_draft.value_counts().unstack(fill_value=0).to_string())
