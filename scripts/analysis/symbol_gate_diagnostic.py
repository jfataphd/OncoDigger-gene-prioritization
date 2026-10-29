"""
symbol_gate_diagnostic.py -- diagnostic of the symbol-presence eligibility filter in run_query_planes.rank_genes (Supplementary Methods S1).

rank_genes reads the engine's gene table in order and skips a gene if it is one of the GENE_EXCLUDE genes, if it was already listed, or if its own
lower-case approved symbol has no posting in the cancer's collection, until 50 genes are listed. This script repeats that step for the primary query
on every cancer with the same engine set-up as run_query_planes.main (frozen root, policy v5, longest-match masking, exact-background hook) and
records, for each engine candidate, whether it was kept or removed and why. It changes nothing and writes one CSV.

Run:    ONCODIGGER_ROOT=<frozen 2026-09-27 root> BGFIX_ENABLE=1 PYTHONHASHSEED=0 python -X utf8 scripts/analysis/symbol_gate_diagnostic.py [out.csv]
Output: data/derived/symbol_gate_diagnostic.csv  (cancer, engine_position, gene, status, final_rank)
        status = kept | GATE_dropped (no posting of its own symbol) | excluded (GENE_EXCLUDE) | duplicate
        The kept rows reproduce data/derived/query_planes/mutation.csv position for position (950 of 950).
"""
import csv
import os
import sys
from pathlib import Path

W = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(W / "scripts" / "analysis"))
import run_query_planes as rq  # noqa: E402

ROOT = Path(os.environ["ONCODIGGER_ROOT"])
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else W / "data" / "derived" / "symbol_gate_diagnostic.csv"
sys.path.insert(0, str(ROOT / "src"))
from oncodigger.search import discovery  # noqa: E402
from oncodigger.search.discovery import associated_entity_tables  # noqa: E402
from oncodigger.search.lexical import bm25_search, normalize_query  # noqa: E402
from oncodigger.search.lexicons import load_lexicons  # noqa: E402

lex = load_lexicons(ROOT)
genes = rq.longest_match_lookups(rq.suppress_ambiguous_symbols(lex["Genes"]), discovery)
lookups = {"Genes": genes}
base = ROOT / "data" / "pubmed_processed"
rows = []
for key, label in rq.CANCER_MAP:
    db = next(p for p in (base / key / "lexical_index.db", base / f"{key}_all" / "lexical_index.db") if p.exists())
    con = rq.open_corpus(db, discovery)
    tokens = normalize_query(rq.PLANES["mutation"])
    hits = bm25_search(con, tokens, rq.SEARCH_POOL, rq.OPTIONS)[:rq.ANALYSIS_PAPERS]
    tables = associated_entity_tables(con, hits, rq.PLANES["mutation"], (2000, 2026), tokens, lookups,
                                      limit=rq.LIMIT, discovery_mode="Primary Evidence", corpus_query="")
    df = tables.get("Genes")
    kept = 0
    seen = set()
    for i, (_, r) in enumerate(df.iterrows(), 1):
        g = rq.GENE_REMAP.get(r["Gene"], r["Gene"])
        status = "kept"
        if g in rq.GENE_EXCLUDE:
            status = "excluded"
        elif g in seen:
            status = "duplicate"
        else:
            n = con.execute("SELECT COUNT(DISTINCT pmid) FROM postings WHERE token = ?", (g.lower(),)).fetchone()[0]
            if n == 0:
                status = "GATE_dropped"
            else:
                seen.add(g)
                kept += 1
        rows.append({"cancer": label, "engine_position": i, "gene": g, "status": status,
                     "final_rank": kept if status == "kept" else ""})
        if kept >= rq.TOP_K and status == "kept":
            break
    con.close()
    d = [r for r in rows if r["cancer"] == label and r["status"] == "GATE_dropped"]
    print(f"{label:20} engine rows read {sum(1 for r in rows if r['cancer'] == label)}  gate-dropped {len(d)}  "
          f"{[(r['gene'], r['engine_position']) for r in d][:8]}", flush=True)
OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["cancer", "engine_position", "gene", "status", "final_rank"])
    w.writeheader()
    w.writerows(rows)
print("wrote", OUT)
