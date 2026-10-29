"""
vocab_corrected_rerun_met_strict.py -- stricter MET robustness condition (Methods 3.3).

The standard correction (vocab_corrected_rerun.py) removes MET mentions that fail the frozen two-step
rule (capitalized standalone match, then no pre-identified alternate meaning). Because "passes both
checks" means only that no counter-evidence was found, a stricter condition keeps MET in Lung Cancer
only where a posting has BOTH a capitalized standalone "MET" and at least one positive gene-context term
(the "genuine" class of met_disambiguation_check.classify), and redacts MET from every other posting.

This script builds that stricter corrected copy (working copy under
data/derived/intogen_rescoring/corrected_corpora_met_strict/, originals untouched), re-runs the real
engine on it exactly as run_query_planes.py does, and records MET's rank and supporting-paper count
together with the resulting Lung Cancer top 10.

Inputs:  ONCODIGGER_ROOT (corpora + raw abstracts)
Output:  data/derived/intogen_rescoring/vocab_corrected_met_strict.csv

Run: python -X utf8 scripts/analysis/vocab_corrected_rerun_met_strict.py
"""
from __future__ import annotations

import csv
import json
import os
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

ONCODIGGER_ROOT = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
STRICT_DIR = ROOT / "data" / "derived" / "intogen_rescoring" / "corrected_corpora_met_strict"
OUT = ROOT / "data" / "derived" / "intogen_rescoring" / "vocab_corrected_met_strict.csv"


def main() -> None:
    if os.environ.get("PYTHONHASHSEED") != "0":
        import subprocess
        sys.exit(subprocess.run([sys.executable, "-X", "utf8", __file__], env=dict(os.environ, PYTHONHASHSEED="0")).returncode)

    import vocab_corrected_rerun as vcr
    from met_disambiguation_check import CAPS_MET, classify
    from run_query_planes import (suppress_ambiguous_symbols, longest_match_lookups, open_corpus, PLANES,
                                  SEARCH_POOL, OPTIONS, ANALYSIS_PAPERS, LIMIT, GENE_REMAP, GENE_EXCLUDE)

    STRICT_DIR.mkdir(parents=True, exist_ok=True)
    vcr.OUT_DIR = STRICT_DIR  # build the working copy elsewhere; the standard corrected copy is untouched

    src_db = ONCODIGGER_ROOT / "data" / "pubmed_processed" / "lung_cancer_inclusive_all" / "lexical_index.db"
    raw = ONCODIGGER_ROOT / "data" / "pubmed_raw" / "lung_cancer_inclusive" / "abstracts.jsonl"
    con = sqlite3.connect(f"file:{src_db.as_posix()}?mode=ro", uri=True)
    met_pmids = {str(p[0]) for p in con.execute("SELECT DISTINCT pmid FROM postings WHERE token = 'met'")}
    con.close()
    def compute_keep():
        keep = set()
        with open(raw, encoding="utf-8") as f:
            for line in f:
                rec = json.loads(line)
                pmid = str(rec.get("pmid", ""))
                if pmid in met_pmids:
                    text = (rec.get("title") or "") + " " + (rec.get("abstract") or "")
                    if CAPS_MET.search(text) and classify(text) == "genuine":
                        keep.add(pmid)
        return keep

    keep = vcr.resolve_pmid_set("met_strict_keep_met_lung_cancer_inclusive", compute_keep, raw)
    print(f"MET postings {len(met_pmids)}; strict condition keeps {len(keep)} ({100 * len(keep) / len(met_pmids):.1f}%)")
    dst = vcr.make_corrected_db("Lung Cancer", "MET", met_pmids - keep)

    sys.path.insert(0, str(ONCODIGGER_ROOT / "src"))
    from oncodigger.search import discovery
    from oncodigger.search.discovery import associated_entity_tables
    from oncodigger.search.lexical import bm25_search, normalize_query
    from oncodigger.search.lexicons import load_lexicons

    lex = load_lexicons(ONCODIGGER_ROOT)
    genes = longest_match_lookups(suppress_ambiguous_symbols(lex["Genes"]), discovery)
    query = PLANES["mutation"]
    c = open_corpus(dst, discovery)
    tokens = normalize_query(query)
    hits = bm25_search(c, tokens, SEARCH_POOL, OPTIONS)[:ANALYSIS_PAPERS]
    df = associated_entity_tables(c, hits, query, (2000, 2026), tokens, {"Genes": genes}, limit=LIMIT,
                                  discovery_mode="Primary Evidence", corpus_query="")["Genes"]
    col = next(x for x in df.columns if x.lower().replace("_", " ") == "supporting papers")
    ranked, seen, met = [], set(), None
    for _, r in df.iterrows():
        g = GENE_REMAP.get(r["Gene"], r["Gene"])
        if g in GENE_EXCLUDE or g in seen:
            continue
        if c.execute("SELECT COUNT(DISTINCT pmid) FROM postings WHERE token = ?", (g.lower(),)).fetchone()[0] == 0:
            continue
        seen.add(g)
        ranked.append(g)
        if g == "MET":
            met = (len(ranked), int(r[col]))
        if len(ranked) >= 50:
            break
    c.close()
    print("Lung top-10 under the strict condition:", ranked[:10])
    print("MET (rank, supporting papers):", met)
    with OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["condition", "met_postings_retained", "met_rank", "met_supporting_papers", "lung_top10"])
        w.writerow(["strict (capitalized + positive context)", len(keep), met[0] if met else "", met[1] if met else 0,
                    ";".join(ranked[:10])])
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
