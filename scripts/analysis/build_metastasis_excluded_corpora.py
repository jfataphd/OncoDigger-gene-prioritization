"""
build_metastasis_excluded_corpora.py -- copies of the collections (all 19, or those named with --corpora) without the records indexed under the
cancer MeSH term's /secondary subheading, with every corpus-wide statistic recomputed on the reduced set (post hoc sensitivity analysis,
Supplementary Results S16).

The removal list is data/derived/metastatic_exclusion/secondary_pmids_<collection>.csv (fetch_secondary_pmids.py), restricted to records that are in
the frozen 2026-09-27 collection. Corpus size N, mean document length, the BM25 IDF table and the postings used for the enrichment background are
rebuilt exactly as build_date_truncated_corpora.py does (same formulas), so run_query_planes.py sees a self-consistent smaller corpus.

Usage: python -X utf8 scripts/analysis/build_metastasis_excluded_corpora.py --src-root <frozen root> --out-root <empty dir>
Output: <out-root>/data/pubmed_processed/<corpus>_cancer_inclusive_all/lexical_index.db and build_report_metastasis_excluded.csv
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import sqlite3
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("bdt", ROOT / "scripts" / "analysis" / "build_date_truncated_corpora.py")
bdt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bdt)
COLLECTIONS = ["breast", "lung", "colorectal", "prostate", "melanoma", "bladder", "kidney", "pancreas", "liver", "stomach", "esophageal",
               "ovarian", "endometrial", "cervical", "thyroid", "brain", "leukemia", "lymphoma", "myeloma"]
LISTS = ROOT / "data" / "derived" / "metastatic_exclusion"


def build_one(name: str, src_root: Path, out_root: Path) -> dict:
    t0 = time.time()
    exclude = [r["pmid"] for r in csv.DictReader((LISTS / f"secondary_pmids_{name}.csv").open(encoding="utf-8")) if r["in_frozen_collection"] == "True"]
    out = out_root / "data" / "pubmed_processed" / f"{name}_cancer_inclusive_all" / "lexical_index.db"
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        out.unlink()
    con = sqlite3.connect("file:" + str(out).replace("\\", "/") + "?mode=rwc", uri=True)
    con.executescript(bdt.SCHEMA)
    con.execute("PRAGMA synchronous=OFF")
    con.execute("PRAGMA journal_mode=OFF")
    con.execute("CREATE TEMP TABLE excl (pmid TEXT PRIMARY KEY)")
    con.executemany("INSERT OR IGNORE INTO excl VALUES (?)", [(p,) for p in exclude])
    con.execute("ATTACH DATABASE ? AS src", ("file:" + str(bdt.src_db(src_root, name)).replace("\\", "/") + "?mode=ro",))
    total = con.execute("SELECT COUNT(*) FROM src.documents").fetchone()[0]
    con.execute("INSERT INTO documents SELECT * FROM src.documents WHERE pmid NOT IN (SELECT pmid FROM excl)")
    n = con.execute("SELECT COUNT(*) FROM documents").fetchone()[0]
    con.execute("INSERT INTO postings SELECT p.token, p.pmid, p.tf FROM src.postings p JOIN documents d ON d.pmid = p.pmid")
    con.execute("CREATE INDEX idx_postings_token ON postings(token)")
    avg = con.execute("SELECT AVG(doc_length) FROM documents").fetchone()[0]
    con.execute("INSERT INTO corpus_stats VALUES (?, ?)", (avg, n))
    rows = [(t, bdt.idf_value(n, df)) for t, df in con.execute("SELECT token, COUNT(*) FROM postings GROUP BY token")]
    con.executemany("INSERT INTO idf VALUES (?, ?)", rows)
    con.commit()
    con.execute("DETACH DATABASE src")
    con.close()
    return dict(corpus=name, source_docs=total, removed_requested=len(exclude), kept_docs=n, removed_actual=total - n,
                idf_tokens=len(rows), avg_doc_length=avg, seconds=round(time.time() - t0))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src-root", type=Path, required=True)
    ap.add_argument("--out-root", type=Path, required=True)
    ap.add_argument("--corpora", nargs="+", default=COLLECTIONS)
    a = ap.parse_args()
    res = []
    for n in a.corpora:
        r = build_one(n, a.src_root, a.out_root)
        print(r, flush=True)
        res.append(r)
    tag = "" if a.corpora == COLLECTIONS else "_" + "_".join(a.corpora)[:40]
    rep = a.out_root / f"build_report_metastasis_excluded{tag}.csv"
    with rep.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(res[0]))
        w.writeheader()
        w.writerows(res)
    (LISTS / rep.name).write_text(rep.read_text(encoding="utf-8"), encoding="utf-8")
    print("wrote", rep)


if __name__ == "__main__":
    main()
