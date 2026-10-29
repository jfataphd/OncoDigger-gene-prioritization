"""
vocab_correction_supporting_papers.py -- supporting-paper counts and ranks for MET (Lung Cancer),
REST (Thyroid Cancer) and FH (Kidney Cancer) before and after the false-positive correction
(Methods 3.3).

Methods 3.3 quotes MET's Lung Cancer supporting-paper count before and after the correction. Those counts
had no recorded generating script. This script runs the real scoring engine, set up exactly as in
run_query_planes.py / vocab_corrected_rescoring.py, on (a) the original corpus database and (b) the
corrected copy built by vocab_corrected_rerun.py (false-positive mentions redacted from title and
abstract text), and records each gene's rank and "Supporting papers" value from the engine's own gene
table.

Inputs:  ONCODIGGER_ROOT (original corpora), data/derived/intogen_rescoring/corrected_corpora/
Output:  data/derived/intogen_rescoring/vocab_correction_supporting_papers.csv

Run: python -X utf8 scripts/analysis/vocab_correction_supporting_papers.py
"""
from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

ONCODIGGER_ROOT = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
CORRECTED_DIR = ROOT / "data" / "derived" / "intogen_rescoring" / "corrected_corpora"
OUT = ROOT / "data" / "derived" / "intogen_rescoring" / "vocab_correction_supporting_papers.csv"

# (corpus key, cancer label, gene)
TARGETS = [
    ("lung_cancer_inclusive", "Lung Cancer", "MET"),
    ("thyroid_cancer_inclusive", "Thyroid Cancer", "REST"),
    ("kidney_cancer_inclusive", "Kidney Cancer", "FH"),
]


def main() -> None:
    if os.environ.get("PYTHONHASHSEED") != "0":
        import subprocess
        env = dict(os.environ, PYTHONHASHSEED="0")
        sys.exit(subprocess.run([sys.executable, "-X", "utf8", __file__], env=env).returncode)

    from run_query_planes import (  # noqa: E402
        suppress_ambiguous_symbols, longest_match_lookups, open_corpus, PLANES, SEARCH_POOL, OPTIONS,
        ANALYSIS_PAPERS, LIMIT, GENE_REMAP, GENE_EXCLUDE,
    )

    sys.path.insert(0, str(ONCODIGGER_ROOT / "src"))
    from oncodigger.search import discovery
    from oncodigger.search.discovery import associated_entity_tables
    from oncodigger.search.lexical import bm25_search, normalize_query
    from oncodigger.search.lexicons import load_lexicons

    lexicons = load_lexicons(ONCODIGGER_ROOT)
    genes = longest_match_lookups(suppress_ambiguous_symbols(lexicons["Genes"]), discovery)
    query = PLANES["mutation"]
    rows = []
    for key, cancer, gene in TARGETS:
        for state, db in (("original", ONCODIGGER_ROOT / "data" / "pubmed_processed" / f"{key}_all" / "lexical_index.db"),
                          ("corrected", CORRECTED_DIR / f"{key}_all" / "lexical_index.db")):
            con = open_corpus(db, discovery)
            tokens = normalize_query(query)
            hits = bm25_search(con, tokens, SEARCH_POOL, OPTIONS)[:ANALYSIS_PAPERS]
            tables = associated_entity_tables(con, hits, query, (2000, 2026), tokens, {"Genes": genes},
                                              limit=LIMIT, discovery_mode="Primary Evidence", corpus_query="")
            df = tables["Genes"]
            rank, seen, found = 0, set(), None
            for _, r in df.iterrows():
                g = GENE_REMAP.get(r["Gene"], r["Gene"])
                if g in GENE_EXCLUDE or g in seen:
                    continue
                try:  # same zero-postings skip as run_query_planes.rank_genes
                    if con.execute("SELECT COUNT(DISTINCT pmid) FROM postings WHERE token = ?", (g.lower(),)).fetchone()[0] == 0:
                        continue
                except Exception:
                    continue
                seen.add(g)
                rank += 1
                if g == gene:
                    found = (rank, r)
                    break
            con.close()
            if found is None:
                rows.append({"cancer": cancer, "gene": gene, "state": state, "rank": "", "supporting_papers": 0})
            else:
                col = next(c for c in df.columns if c.lower().replace("_", " ") == "supporting papers")
                rows.append({"cancer": cancer, "gene": gene, "state": state, "rank": found[0],
                             "supporting_papers": int(found[1][col])})
            print(rows[-1])
    with OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["cancer", "gene", "state", "rank", "supporting_papers"])
        w.writeheader()
        w.writerows(rows)
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
