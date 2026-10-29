"""
extract_features.py -- per-gene components of the OncoDigger relevance score, for formula-sensitivity work (exploratory).

For each of the 19 cancers, runs the paper's gene-ranking path (same query, pool, vocabulary treatment and options as
scripts/analysis/run_query_planes.py) but asks for every gene with at least 2 supporting papers instead of the top 100, and records
the quantities the relevance formula is built from:

    S  supporting papers in the pool            B  background papers used by the engine (its 'Corpus papers' column)
    N  pool size (retrieved papers)             C  papers in the cancer collection
    W  rank-weighted evidence = sum over supporting papers of 1/sqrt(r), r = 1-based pool rank (the engine uses 1/sqrt(zero-based rank + 1))
    ranks  the pool ranks (1-based) of the supporting papers, for alternative rank decays

W is recomputed here from the supporting PMIDs and the pool order and checked against the engine's displayed Relevance
(relevance = W * ln(1+E) * ln(1+S), E = (S/N)/(B/C)). The baseline top 50 is checked against data/derived/query_planes/mutation_vocab_corrected.csv.

Read-only. Run: PYTHONHASHSEED=0 python -X utf8 extract_features.py <ONCODIGGER_ROOT>   (use the corrected root, scratchpad s27corr)
Output: features_raw.csv (cancer, gene, baseline_rank, S, B, N, C, W, ranks, displayed_relevance)
"""
import importlib.util
import math
import os
import sys
from pathlib import Path

import pandas as pd

assert os.environ.get("PYTHONHASHSEED") == "0", "set PYTHONHASHSEED=0"
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ROOT = Path(sys.argv[1])
sys.path.insert(0, str(REPO / "scripts" / "analysis"))
spec = importlib.util.spec_from_file_location("rqp", REPO / "scripts" / "analysis" / "run_query_planes.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

sys.path.insert(0, str(ROOT / "src"))
from oncodigger.search import discovery
from oncodigger.search.discovery import associated_entity_tables
from oncodigger.search.lexical import bm25_search, normalize_query
from oncodigger.search.lexicons import load_lexicons

lex = load_lexicons(ROOT)
genes = m.longest_match_lookups(m.suppress_ambiguous_symbols(lex["Genes"]), discovery)
QUERY = m.PLANES["mutation"]
tokens = normalize_query(QUERY)

rows = []
for key, label in m.CANCER_MAP:
    db = ROOT / "data" / "pubmed_processed" / f"{key}_all" / "lexical_index.db"
    con = m.open_corpus(db, discovery)
    hits = bm25_search(con, tokens, m.SEARCH_POOL, m.OPTIONS)[: m.ANALYSIS_PAPERS]
    pool_rank = {str(h["pmid"]): i + 1 for i, h in enumerate(hits)}
    N = len(hits)
    C = con.execute("SELECT COUNT(*) FROM documents").fetchone()[0]
    tabs = associated_entity_tables(con, hits, QUERY, (2000, 2026), tokens, {"Genes": genes}, limit=100000,
                                    discovery_mode="Primary Evidence", corpus_query="")
    g = tabs["Genes"]
    seen, k = set(), 0
    for _, r in g.iterrows():
        gene = m.GENE_REMAP.get(r["Gene"], r["Gene"])
        if gene in m.GENE_EXCLUDE or gene in seen:
            continue
        try:
            n = con.execute("SELECT COUNT(DISTINCT pmid) FROM postings WHERE token = ?", (gene.lower(),)).fetchone()[0]
        except Exception:
            n = 0
        if n == 0:
            continue
        seen.add(gene)
        k += 1
        pm = [p for p in str(r["Supporting PMIDs"]).split(";") if p]
        rk = sorted(pool_rank[p] for p in pm if p in pool_rank)
        rows.append({"cancer": label, "gene": gene, "baseline_rank": k, "S": int(r["Supporting papers"]), "n_listed": len(rk),
                     "B": int(r["Corpus papers"]), "N": N, "C": C, "W": sum(1 / math.sqrt(x) for x in rk),
                     "ranks": ";".join(map(str, rk)), "displayed_relevance": float(r["Relevance"])})
    con.close()
    print(label, N, C, k, flush=True)
pd.DataFrame(rows).to_csv(HERE / "features_raw.csv", index=False)
print("wrote features_raw.csv", len(rows))
