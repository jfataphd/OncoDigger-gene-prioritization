"""Read-only extraction of the pool abstracts crediting each top-10 gene, any cancer. usage: extract_any.py "<Cancer label>" <key>
Frozen engine + frozen v5 policy; supporting pool papers do not depend on the background count."""
import os, sys
from pathlib import Path
os.environ.setdefault("PYTHONHASHSEED", "0")
label, key = sys.argv[1], sys.argv[2]
REPO = Path(os.environ.get("ONCODIGGER_PAPER_REPO") or Path(__file__).resolve().parents[3]); HERE = Path(__file__).resolve().parent; ROOT = HERE/"root"
sys.path.insert(0, str(REPO/"scripts/analysis")); sys.path.insert(0, str(ROOT/"src"))
import run_query_planes as rq
from oncodigger.search import discovery
from oncodigger.search.lexical import bm25_search, normalize_query
from oncodigger.search.lexicons import load_lexicons
import pandas as pd
TOP10 = pd.read_csv(REPO/"data/canonical/precision_at_k_detail.csv").query("cancer==@label and rank<=10").gene.tolist()
genes = rq.longest_match_lookups(rq.suppress_ambiguous_symbols(load_lexicons(ROOT)["Genes"]), discovery)
db = ROOT/f"data/pubmed_processed/{key}_cancer_inclusive_all/lexical_index.db"
con = rq.open_corpus(db, discovery)
q = rq.PLANES["mutation"]; toks = normalize_query(q)
hits = bm25_search(con, toks, rq.SEARCH_POOL, rq.OPTIONS)[:rq.ANALYSIS_PAPERS]
tables = discovery.associated_entity_tables(con, hits, q, (2000,2026), toks, {"Genes": genes}, limit=100, discovery_mode="Primary Evidence", corpus_query="")
g = tables["Genes"]; query_terms = {t.replace("-"," ") for t in toks}; docs = discovery._docs_for_results(con, hits)
rows=[]; summary={}
for gene in TOP10:
    r = g[g["Gene"].isin([gene, {"TP53":"P53","ERBB2":"HER2"}.get(gene,gene)])]
    pm=set()
    for _,x in r.iterrows(): pm |= set(x["Supporting PMIDs"].split(";"))
    summary[gene]=len(pm)
    singles=[t for t,es in genes["single"].items() if any(e["label"]==gene for e in es)]
    phrases=[p for a,l in genes["phrases"].items() for p,e in l if e["label"]==gene]
    for rank,h in enumerate(hits,1):
        if h["pmid"] not in pm: continue
        d = docs[h["pmid"]]; text=discovery.full_document_text(d); padded=f" {text} "; tset=set(discovery.document_tokens(d))
        m=[t for t in singles if t in tset and t not in query_terms]+[p for p in phrases if p not in query_terms and f" {p} " in padded]
        rows.append(dict(gene=gene,pmid=h["pmid"],pool_rank=rank,year=d.get("year") or d.get("publication_year"),matched_terms="|".join(sorted(m)),title=d.get("title",""),abstract=d.get("abstract","")))
out=HERE/f"{key}_pool_gene_matches.csv"; pd.DataFrame(rows).to_csv(out,index=False,encoding="utf-8")
print(label,"pool",len(hits),summary,"rows",len(rows))
