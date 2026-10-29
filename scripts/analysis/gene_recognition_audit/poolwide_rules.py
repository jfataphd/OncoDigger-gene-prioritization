"""Whole-pool automated counts for each recurring error mechanism (text rules only; each rule was validated on abstracts the reader read in full).
Counts are of abstracts credited to the gene that satisfy the mechanism's rule. Compared with the number the reader judged wrong in the read set."""
import pandas as pd, re
K={"breast":"Breast Cancer","lung":"Lung Cancer","colorectal":"Colorectal Cancer","prostate":"Prostate Cancer","melanoma":"Melanoma","bladder":"Bladder Cancer","kidney":"Kidney Cancer","pancreas":"Pancreatic Cancer","liver":"Liver Cancer","stomach":"Stomach Cancer","esophageal":"Esophageal Cancer","ovarian":"Ovarian Cancer","endometrial":"Endometrial Cancer","cervical":"Cervical Cancer","thyroid":"Thyroid Cancer","brain":"Brain Cancer","leukemia":"Leukemia","lymphoma":"Lymphoma","myeloma":"Myeloma"}
VAR={"MLH1":r"h?mlh[ -]?1|mutl homolog 1","MSH2":r"h?msh[ -]?2|muts homolog 2"}
BCELL=re.compile(r"b[- ]cell (antigen )?receptor|bcr signal|bcr[- ]stereotyp|bcr[- ]?(pathway|inhibitor)|surface bcr"); ABL=re.compile(r"bcr[- /:]{1,2}abl|breakpoint cluster|philadelphia|\bph\+|abl1")
CHEM=re.compile(r"tert[- ]?(?:butyl|butoxy|amyl|octyl|pentyl)|tertiary butyl"); APCC=re.compile(r"apc[ /-]c\b|anaphase[- ]promoting|antigen[- ]presenting|activated protein c|argon plasma|\bapcs\b")
rows=[]
for k,l in K.items():
    d=pd.read_csv(f"{k}_pool_gene_matches.csv",dtype=str).fillna("")
    for g,x in d.groupby("gene"):
        t=(x.title+" "+x.abstract).str.replace("_"," ").str.lower(); terms=x.matched_terms.str.split("|"); sym=g.lower(); n=len(x)
        def add(mech,cnt): 
            if cnt: rows.append(dict(cancer=l,gene=g,mechanism=mech,pool_credits=n,n_pool_rule=int(cnt),pct=round(100*cnt/n,1)))
        if g=="PIK3CA": add("PI3K -> PIK3CA",sum(1 for tt in terms if tt==["pi3k"]))
        if g in VAR: add("HNPCC syndrome name",sum(1 for tt,te in zip(terms,t) if "hnpcc" in tt and sym not in tt and not re.search(r"(?<![a-z0-9])(%s)(?![a-z0-9])"%VAR[g],te)))
        if g=="EGFR": add("EGFR via HER2 name",sum(1 for tt,te in zip(terms,t) if sym not in tt and re.search(r"receptor[ -]?2\b",te)))
        if g=="SETD2": add("HIF-1 -> SETD2",sum(1 for tt in terms if tt==["hif 1"]))
        if g=="WT1": add("WAGR -> WT1",sum(1 for tt in terms if tt==["wagr"]))
        if g=="PARP1": add("PARP family/drug-class -> PARP1",sum(1 for tt,te in zip(terms,t) if tt==["parp"] and not re.search(r"parp[ -]?1\b",te)))
        if g=="TERT": add("tert-butyl chemistry",sum(1 for te in t if CHEM.search(te) and not re.search(r"(?<![a-z0-9])tert(?![a-z0-9])|telomerase reverse transcriptase",CHEM.sub(" ",te))))
        if g=="APC": add("APC/C or other non-gene APC",sum(1 for te in t if APCC.search(te) and not re.search(r"adenomatous|polyposis|\bfap\b|apc (mutation|gene|loss|variant|allele|germline|somatic|alteration|status)|wnt|catenin|apc [a-z]\.\d|apc c ?\d",te)))
        if g=="BCR": add("BCR = B-cell receptor",sum(1 for te in t if BCELL.search(te) and not ABL.search(te)))
        if g=="CD274": add("CD274 credited only via 'PD-L1' (symbol absent)",sum(1 for tt in terms if sym not in tt))
        if g=="ERBB2": add("ERBB2 credited only via HER2/HER-2 (symbol absent)",sum(1 for tt in terms if "erbb2" not in tt and "erb_b2_receptor_tyrosine_kinase_2" not in tt))
        if g=="CTNNB1": add("CTNNB1 credited only via 'beta-catenin' (symbol absent; same protein/gene)",sum(1 for tt in terms if sym not in tt))
R=pd.DataFrame(rows); R.to_csv("ALL19_poolwide_rule_counts.csv",index=False); print(R.sort_values("pct",ascending=False).to_string(index=False))
