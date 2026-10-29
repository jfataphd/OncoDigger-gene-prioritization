"""Read-only Breast counterfactual: exact-B ranking of the Breast pool under alternative gene-term sets. Writes only to this folder."""
import os, sys, json
from pathlib import Path
REPO = Path(os.environ.get("ONCODIGGER_PAPER_REPO") or Path(__file__).resolve().parents[3])
HERE = Path(__file__).resolve().parent; ROOT = HERE/"root"
os.environ.update(PYTHONHASHSEED="0", PYTHONUTF8="1", BGFIX_ENABLE="1", BGFIX_CACHE=str(HERE/"bgfix_cache"), BGFIX_WORKERS="8",
                  ONCODIGGER_ROOT=str(ROOT), GENE_QUALITY_TEST_ENGINE=str(REPO/"reference_engine_code/gene_quality_test_engine"))
(HERE/"bgfix_cache").mkdir(exist_ok=True)
if "bgfix_hook" not in sys.modules and os.environ.get("_CF_CHILD")!="1":
    import subprocess
    env=dict(os.environ, _CF_CHILD="1", PYTHONPATH=str(REPO/"exploratory/background_rebuild/hook"))
    sys.exit(subprocess.run([sys.executable,"-X","utf8"]+sys.argv,env=env).returncode)
sys.path.insert(0,str(REPO/"scripts/analysis")); sys.path.insert(0,str(ROOT/"src"))
import run_query_planes as rq, copy
from oncodigger.search import discovery
from oncodigger.search.lexical import bm25_search, normalize_query
from oncodigger.search.lexicons import load_lexicons
variant=sys.argv[1]
lex=load_lexicons(ROOT)
genes=rq.longest_match_lookups(rq.suppress_ambiguous_symbols(lex["Genes"]), discovery)
genes={"single":{k:list(v) for k,v in genes["single"].items()},"phrases":{k:list(v) for k,v in genes["phrases"].items()}}
def drop(labels, singles=(), phrases=()):
    for t in singles:
        if t in genes["single"]:
            genes["single"][t]=[e for e in genes["single"][t] if e["label"] not in labels]
            if not genes["single"][t]: del genes["single"][t]
    for a,l in list(genes["phrases"].items()):
        genes["phrases"][a]=[(p,e) for p,e in l if not (p in phrases and e["label"] in labels)]
if variant in ("erbb2_strict","erbb2_strict_no_pi3k"): drop({"ERBB2","HER2"}, singles=("her2",), phrases=("her 2",))
if variant=="erbb2_strict_no_pi3k": drop({"PIK3CA"}, singles=("pi3k",), phrases=("pi3k",))
for lab in ("ERBB2","HER2","PIK3CA"):
    print(lab,"terms:",sorted({t for t,es in genes["single"].items() if any(e["label"]==lab for e in es)}|{p for a,l in genes["phrases"].items() for p,e in l if e["label"]==lab}),flush=True)
con=rq.open_corpus(ROOT/"data/pubmed_processed/breast_cancer_inclusive_all/lexical_index.db", discovery)
q=rq.PLANES["mutation"]; toks=normalize_query(q)
hits=bm25_search(con,toks,rq.SEARCH_POOL,rq.OPTIONS)[:rq.ANALYSIS_PAPERS]
t=discovery.associated_entity_tables(con,hits,q,(2000,2026),toks,{"Genes":genes},limit=rq.LIMIT,discovery_mode="Primary Evidence",corpus_query="")["Genes"]
out=[];seen=set()
for _,r in t.iterrows():
    g=rq.GENE_REMAP.get(r["Gene"],r["Gene"])
    if g in rq.GENE_EXCLUDE or g in seen: continue
    n=con.execute("SELECT COUNT(DISTINCT pmid) FROM postings WHERE token=?",(g.lower(),)).fetchone()[0]
    if n==0: continue
    seen.add(g); out.append(dict(rank=len(out)+1,gene=g,relevance=float(r["Relevance"]),support=int(r["Supporting papers"]),background=int(r["Corpus papers"])))
    if len(out)>=50: break
json.dump(out,open(HERE/f"cf_{variant}.json","w"),indent=1)
print(variant,[o["gene"] for o in out[:15]],flush=True)
