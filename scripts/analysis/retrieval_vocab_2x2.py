"""
retrieval_vocab_2x2.py -- manuscript1.1 (IntOGen-primary draft) ONLY, diagnostic/construction
work, not yet wired into the manuscript. Decomposes the ~48-point gap between the matched
frequency-only baseline (BM25 1,000-paper pool + fully cleaned vocabulary + flat count =
91.6% any-cancer / 76.3% exact concordance, scripts/analysis/matched_frequency_baseline.py)
and the old full-corpus baseline (43.2%) into its two confounded factors: BM25 retrieval
(1,000-paper pool vs. full corpus) and vocabulary cleanliness (less-cleaned vs. fully
cleaned), via a 2x2 factorial design, all four cells computed fresh and consistently.

Method: single-word gene-symbol/alias token matching directly against each corpus's SQLite
`postings` inverted index (COUNT DISTINCT pmid per gene, unioning across all of a gene's
single-word tokens), NOT the full associated_entity_tables() Python entity-matching pipeline.
This is a documented simplification: it omits multi-word phrase aliases and the longest-match
nested-name correction, which the real pipeline's enrichment/relevance path includes. It is
adequate for this decomposition because (a) the 91.6%/76.3% matched-baseline cell this design
reproduces as a sanity check should come close to (not necessarily identical to) the real
associated_entity_tables()-based result, and (b) the GENE_EXCLUDE/suppress_ambiguous_symbols
distinction (the actual "less-cleaned vs fully cleaned" axis under test) acts on single-word
symbol tokens specifically, so this simplification does not remove the signal being measured.

Four cells:
  Cell 1: full corpus (every pmid in the corpus), less-cleaned vocabulary (raw load_lexicons()
          single-token map, no GENE_EXCLUDE, no suppress_ambiguous_symbols).
  Cell 2: full corpus, fully cleaned vocabulary (GENE_EXCLUDE + suppress_ambiguous_symbols
          applied to the single-token map).
  Cell 3: BM25 1,000-paper pool (same retrieval as the real pipeline: mutation/mutations
          query, SEARCH_POOL=1250 candidates, ANALYSIS_PAPERS=1000), less-cleaned vocabulary.
  Cell 4: BM25 1,000-paper pool, fully cleaned vocabulary (closest cell to the existing
          91.6%/76.3% matched_frequency_baseline.py result; reproduced here fresh for a
          consistent four-way comparison, NOT by reusing that script's output file).

Tie-breaking: alphabetical gene symbol (independent of the real relevance formula), fixing
the tie-breaking leak in the earlier matched_frequency_baseline.py (which broke ties using
the real formula's own relevance order, a stable-sort artifact that left the "frequency-only"
result not fully independent of the formula it was meant to isolate).

Outputs (new files only; does not touch or overwrite matched_frequency_baseline_detail.csv,
matched_frequency_baseline_summary.csv, or freq_baseline_top50_aliases.csv):
  data/derived/intogen_rescoring/retrieval_vocab_2x2_detail.csv
  data/derived/intogen_rescoring/retrieval_vocab_2x2_summary.csv
  data/derived/intogen_rescoring/retrieval_vocab_2x2_cell4_tie_audit.csv

Run: python -X utf8 scripts/analysis/retrieval_vocab_2x2.py [--cancers "Lung Cancer" ...]
"""
from __future__ import annotations

import argparse
import csv
import os
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ONCODIGGER = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
OUT_DIR = ROOT / "data" / "derived" / "intogen_rescoring"
INTOGEN_TSV = ROOT / "data" / "external" / "intogen" / "drivers" / "2024-06-18_IntOGen-Drivers" / "Compendium_Cancer_Genes.tsv"

TOP_K = 10
SEARCH_POOL, ANALYSIS_PAPERS = 1250, 1000
MIN_SUPPORT = 2
GENE_EXCLUDE = {"MICE", "DCR", "PC", "C2", "FBN1", "SGCG", "WDHD1"}
SUPPRESS_SYMBOLS = {"GC", "HR", "HCCS", "BPIFA4P", "MB", "ELN", "CP", "UBC", "CAMP", "EFS", "IVD", "SCT", "SON"}
GENE_REMAP = {"P53": "TP53", "HER2": "ERBB2"}
QUERY_TOKENS = ["mutation", "mutations"]
OPTIONS = {
    "year_min": 2000, "year_max": 2026, "exclude_reviews": False, "downrank_reviews": True,
    "title_match_multiplier": 1.2, "mesh_match_multiplier": 1.15, "recency_boost": True,
    "discovery_mode": "Primary Evidence",
}
CANCER_MAP = [
    ("breast_cancer_inclusive", "Breast Cancer"), ("lung_cancer_inclusive", "Lung Cancer"),
    ("colorectal_cancer_inclusive", "Colorectal Cancer"), ("prostate_cancer_inclusive", "Prostate Cancer"),
    ("melanoma_cancer_inclusive", "Melanoma"), ("bladder_cancer_inclusive", "Bladder Cancer"),
    ("kidney_cancer_inclusive", "Kidney Cancer"), ("pancreas_cancer_inclusive", "Pancreatic Cancer"),
    ("liver_cancer_inclusive", "Liver Cancer"), ("stomach_cancer_inclusive", "Stomach Cancer"),
    ("esophageal_cancer_inclusive", "Esophageal Cancer"), ("ovarian_cancer_inclusive", "Ovarian Cancer"),
    ("endometrial_cancer_inclusive", "Endometrial Cancer"), ("cervical_cancer_inclusive", "Cervical Cancer"),
    ("thyroid_cancer_inclusive", "Thyroid Cancer"), ("brain_cancer_inclusive", "Brain Cancer"),
    ("leukemia_cancer_inclusive", "Leukemia"), ("lymphoma_cancer_inclusive", "Lymphoma"),
    ("myeloma_cancer_inclusive", "Myeloma"),
]
INTOGEN_CANCER_MAP = {
    "Breast Cancer": {"BRCA"}, "Lung Cancer": {"LUAD", "LUSC", "NSCLC", "SCLC"},
    "Colorectal Cancer": {"COAD", "READ", "COADREAD"}, "Prostate Cancer": {"PRAD"},
    "Melanoma": {"MEL", "SKCM", "UM"}, "Bladder Cancer": {"BLCA", "UTUC"},
    "Kidney Cancer": {"CCRCC", "CHRCC", "PRCC", "RCC", "WT"}, "Pancreatic Cancer": {"PAAD", "PANET"},
    "Liver Cancer": {"HCC", "CHOL", "LIHB"}, "Stomach Cancer": {"STAD"},
    "Esophageal Cancer": {"ESCA", "ESCC"}, "Ovarian Cancer": {"OVT"},
    "Endometrial Cancer": {"UCEC", "UCS"}, "Cervical Cancer": {"CESC", "CEAD"},
    "Thyroid Cancer": {"WDTC"}, "Brain Cancer": {"GB", "GBM", "HGGNOS", "LGGNOS", "PAST", "MBL", "EPM", "ATRT"},
    "Leukemia": {"ALL", "AML", "CLLSLL", "MDS"}, "Lymphoma": {"BL", "DLBCLNOS", "NHL", "MLYM"},
    "Myeloma": {"PCM"},
}


def build_token_to_gene(single_lexicon: dict, cleaned: bool) -> dict[str, str]:
    """token -> gene symbol, one winner per token (first entity id), optionally with
    GENE_EXCLUDE genes dropped entirely and SUPPRESS_SYMBOLS bare tokens removed."""
    token_to_gene = {}
    for token, entities in single_lexicon.items():
        if not entities:
            continue
        gene = entities[0]["id"]
        gene = GENE_REMAP.get(gene, gene)
        if cleaned and gene in GENE_EXCLUDE:
            continue
        if cleaned and gene in SUPPRESS_SYMBOLS and token == gene.lower():
            continue
        token_to_gene[token] = gene
    return token_to_gene


def bm25_pool_pmids(con: sqlite3.Connection) -> set[str]:
    """Reproduces the real pipeline's BM25 pool selection (k1=1.5, b=0.75, title/MeSH/
    recency multipliers, year filter), identical formula to lexical.py's _bm25_sqlite,
    returning the top ANALYSIS_PAPERS pmids by score."""
    import math
    row = con.execute("SELECT avg_doc_length FROM corpus_stats").fetchone()
    avgdl = row[0] or 1
    k1, b = 1.5, 0.75
    placeholders = ",".join("?" * len(QUERY_TOKENS))
    idf_map = dict(con.execute(
        f"SELECT token, idf FROM idf WHERE token IN ({placeholders})", QUERY_TOKENS).fetchall())

    # Single batched fetch of (pmid -> {token: tf}) for both query tokens, no per-pmid queries.
    pmid_tf: dict[str, dict[str, int]] = defaultdict(dict)
    for token in QUERY_TOKENS:
        if token not in idf_map:
            continue
        for pmid, tf in con.execute("SELECT pmid, tf FROM postings WHERE token = ?", (token,)).fetchall():
            pmid_tf[pmid][token] = tf

    matched_pmids = list(pmid_tf.keys())
    docs = {}
    CHUNK = 500
    for i in range(0, len(matched_pmids), CHUNK):
        chunk = matched_pmids[i:i + CHUNK]
        ph = ",".join("?" * len(chunk))
        rows = con.execute(
            f"SELECT pmid, title, mesh_terms, publication_year, publication_types, doc_length "
            f"FROM documents WHERE pmid IN ({ph})", chunk).fetchall()
        for pmid, title, mesh, year, ptypes, doclen in rows:
            docs[pmid] = (title or "", mesh or "", year, ptypes or "", doclen or avgdl)

    final_scores = {}
    for pmid, tf_by_token in pmid_tf.items():
        if pmid not in docs:
            continue
        title, mesh, year, ptypes, doclen = docs[pmid]
        if OPTIONS["year_min"] and year and year < OPTIONS["year_min"]:
            continue
        if OPTIONS["year_max"] and year and year > OPTIONS["year_max"]:
            continue
        tf_sum = 0.0
        for token, tf in tf_by_token.items():
            idf = idf_map.get(token)
            if idf is None or tf <= 0:
                continue
            denom = tf + k1 * (1 - b + b * doclen / avgdl)
            tf_sum += idf * (tf * (k1 + 1)) / denom
        mult = 1.0
        low_title = title.lower()
        if any(t in low_title for t in QUERY_TOKENS):
            mult *= OPTIONS["title_match_multiplier"]
        if any(t in (mesh or "").lower() for t in QUERY_TOKENS):
            mult *= OPTIONS["mesh_match_multiplier"]
        if OPTIONS["recency_boost"] and year:
            if year >= 2020:
                mult *= 1.1
            elif year < 2010:
                mult *= 0.95
        if "Review" in (ptypes or ""):
            mult *= 0.65
        final_scores[pmid] = tf_sum * mult

    ranked = sorted(final_scores.items(), key=lambda kv: kv[1], reverse=True)[:ANALYSIS_PAPERS]
    return {pmid for pmid, _ in ranked}


def gene_counts(con: sqlite3.Connection, token_to_gene: dict[str, str], pmid_filter: set[str] | None) -> dict[str, int]:
    """Distinct-pmid count per gene, unioning all of that gene's single-word tokens,
    optionally restricted to pmid_filter (the BM25 pool) or None for the full corpus."""
    gene_pmids: dict[str, set[str]] = defaultdict(set)
    tokens = list(token_to_gene.keys())
    CHUNK = 500
    for i in range(0, len(tokens), CHUNK):
        chunk = tokens[i:i + CHUNK]
        ph = ",".join("?" * len(chunk))
        rows = con.execute(f"SELECT token, pmid FROM postings WHERE token IN ({ph})", chunk).fetchall()
        for token, pmid in rows:
            if pmid_filter is not None and pmid not in pmid_filter:
                continue
            gene = token_to_gene.get(token)
            if gene:
                gene_pmids[gene].add(pmid)
    return {g: len(pmids) for g, pmids in gene_pmids.items() if len(pmids) >= MIN_SUPPORT}


def top10_alpha_tiebreak(counts: dict[str, int]) -> list[str]:
    ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    return [g for g, _ in ranked[:TOP_K]]


def load_intogen():
    genes_any, pairs = set(), set()
    with open(INTOGEN_TSV, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sym = row["SYMBOL"].strip().upper()
            genes_any.add(sym)
            code = row["CANCER_TYPE"].strip()
            for cancer, codes in INTOGEN_CANCER_MAP.items():
                if code in codes:
                    pairs.add((sym, cancer))
    return genes_any, pairs


def score(top10_by_cancer: dict[str, list[str]], genes_any: set, exact_pairs: set) -> tuple[int, int, int]:
    any_hits = exact_hits = n = 0
    for cancer, genes in top10_by_cancer.items():
        for g in genes:
            n += 1
            gu = g.strip().upper()
            if gu in genes_any:
                any_hits += 1
            if (gu, cancer) in exact_pairs:
                exact_hits += 1
    return any_hits, exact_hits, n


def main() -> None:
    if os.environ.get("PYTHONHASHSEED") != "0":
        import subprocess
        env = dict(os.environ, PYTHONHASHSEED="0")
        sys.exit(subprocess.run([sys.executable, "-X", "utf8"] + sys.argv, env=env).returncode)

    ap = argparse.ArgumentParser()
    ap.add_argument("--oncodigger-root", type=Path, default=DEFAULT_ONCODIGGER)
    ap.add_argument("--cancers", nargs="+", default=None)
    args = ap.parse_args()

    sys.path.insert(0, str(args.oncodigger_root / "src"))
    from oncodigger.search.lexicons import load_lexicons

    lexicons = load_lexicons(args.oncodigger_root)
    single = lexicons["Genes"]["single"]
    token_to_gene_dirty = build_token_to_gene(single, cleaned=False)
    token_to_gene_clean = build_token_to_gene(single, cleaned=True)
    print(f"vocab sizes: less-cleaned={len(token_to_gene_dirty)} tokens, "
          f"fully cleaned={len(token_to_gene_clean)} tokens", flush=True)

    db_base = args.oncodigger_root / "data" / "pubmed_processed"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    cells = {"cell1_fullcorpus_dirty": {}, "cell2_fullcorpus_clean": {},
             "cell3_pool_dirty": {}, "cell4_pool_clean": {}}
    detail_rows = []
    tie_audit_rows = []

    for key, label in CANCER_MAP:
        if args.cancers and label not in args.cancers:
            continue
        db = next((p for p in (db_base / key / "lexical_index.db", db_base / f"{key}_all" / "lexical_index.db")
                   if p.exists()), None)
        if db is None:
            print(f"  MISSING index for {label}, skipping", flush=True)
            continue
        con = sqlite3.connect(str(db))
        print(f"  {label}: building BM25 pool...", flush=True)
        pool = bm25_pool_pmids(con)
        print(f"  {label}: pool size = {len(pool)}", flush=True)

        counts_c1 = gene_counts(con, token_to_gene_dirty, pmid_filter=None)
        counts_c2 = gene_counts(con, token_to_gene_clean, pmid_filter=None)
        counts_c3 = gene_counts(con, token_to_gene_dirty, pmid_filter=pool)
        counts_c4 = gene_counts(con, token_to_gene_clean, pmid_filter=pool)
        con.close()

        for cell_name, counts in [("cell1_fullcorpus_dirty", counts_c1), ("cell2_fullcorpus_clean", counts_c2),
                                   ("cell3_pool_dirty", counts_c3), ("cell4_pool_clean", counts_c4)]:
            top10 = top10_alpha_tiebreak(counts)
            cells[cell_name][label] = top10
            for rank, gene in enumerate(top10, 1):
                detail_rows.append({"cell": cell_name, "cancer": label, "rank": rank, "gene": gene,
                                     "count": counts[gene]})

        # Tie audit for cell 4 specifically: any count tied at/crossing the rank-10 boundary?
        ranked_c4 = sorted(counts_c4.items(), key=lambda kv: (-kv[1], kv[0]))
        if len(ranked_c4) > TOP_K:
            boundary_count = ranked_c4[TOP_K - 1][1]
            tied = [g for g, c in ranked_c4 if c == boundary_count]
            if len(tied) > 1:
                tie_audit_rows.append({"cancer": label, "boundary_count": boundary_count,
                                        "tied_genes": ";".join(sorted(tied)), "n_tied": len(tied)})

        print(f"  {label:20} c1(dirty,full)={' '.join(cells['cell1_fullcorpus_dirty'][label])}", flush=True)
        print(f"  {label:20} c2(clean,full)={' '.join(cells['cell2_fullcorpus_clean'][label])}", flush=True)
        print(f"  {label:20} c3(dirty,pool)={' '.join(cells['cell3_pool_dirty'][label])}", flush=True)
        print(f"  {label:20} c4(clean,pool)={' '.join(cells['cell4_pool_clean'][label])}", flush=True)

    with (OUT_DIR / "retrieval_vocab_2x2_detail.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["cell", "cancer", "rank", "gene", "count"])
        w.writeheader()
        w.writerows(detail_rows)

    with (OUT_DIR / "retrieval_vocab_2x2_cell4_tie_audit.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["cancer", "boundary_count", "tied_genes", "n_tied"])
        w.writeheader()
        w.writerows(tie_audit_rows)

    genes_any, exact_pairs = load_intogen()
    summary = []
    for cell_name in cells:
        any_hits, exact_hits, n = score(cells[cell_name], genes_any, exact_pairs)
        summary.append({"cell": cell_name, "any_cancer_precision_at_10": f"{any_hits}/{n}",
                         "any_cancer_pct": round(100 * any_hits / n, 1) if n else None,
                         "exact_concordance": f"{exact_hits}/{n}",
                         "exact_pct": round(100 * exact_hits / n, 1) if n else None})
        print(f"{cell_name}: any-cancer {any_hits}/{n} = {100*any_hits/n:.1f}%, "
              f"exact {exact_hits}/{n} = {100*exact_hits/n:.1f}%")

    with (OUT_DIR / "retrieval_vocab_2x2_summary.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["cell", "any_cancer_precision_at_10", "any_cancer_pct",
                                           "exact_concordance", "exact_pct"])
        w.writeheader()
        w.writerows(summary)

    print()
    print(f"Tie audit: {len(tie_audit_rows)} of 19 cancers had a boundary tie in cell 4.")


if __name__ == "__main__":
    main()
