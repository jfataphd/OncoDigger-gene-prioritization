"""
build_date_truncated_corpora.py -- build publication-year-truncated copies of the 19 OncoDigger
lexical-index SQLite corpora, with every corpus-wide statistic recomputed on the truncated set.

Why: the engine's year_max option only restricts which papers BM25 retrieves. The BM25 IDF table,
the average document length, the corpus size N and the per-gene background counts used by the
enrichment score are all computed on the WHOLE database, so a year_max run is not a clean
"corpus as of year X" run. This script builds an actual truncated corpus so the full pipeline
(run_query_planes.py, unchanged except for --oncodigger-root) sees only papers up to the cutoff.

Formulas (identical to scripts/build_lexical_index.py and migrate_to_sqlite.py of the OncoDigger
application at commit 4e39fc0):
  doc_count       = number of documents
  avg_doc_length  = mean of documents.doc_length
  df(token)       = number of postings rows for the token
  idf(token)      = ln(1 + (N - df + 0.5) / (df + 0.5)), tokens with df > 0 only
  enrichment background = COUNT(*) of postings rows for the token (read at query time)

Documents with a NULL publication_year are excluded by the engine's own year filter
(`d.publication_year <= ?`), so they are excluded here too, and counted in the report.

--verify mode: recompute idf / corpus_stats from the FULL database and compare them to the stored
tables. A rebuild is only trusted if this agrees with the stored values.

Usage:
  python -X utf8 scripts/analysis/build_date_truncated_corpora.py --src-root <OncoDigger root> \
      --out-root <empty dir> --year-max 2022
  python -X utf8 scripts/analysis/build_date_truncated_corpora.py --src-root <OncoDigger root> --verify
Output: <out-root>/data/pubmed_processed/<corpus>/lexical_index.db and build_report.csv
"""
from __future__ import annotations

import argparse
import csv
import math
import sqlite3
import sys
import time
from pathlib import Path

CORPORA = [
    "breast", "lung", "colorectal", "prostate", "melanoma", "bladder", "kidney", "pancreas", "liver",
    "stomach", "esophageal", "ovarian", "endometrial", "cervical", "thyroid", "brain", "leukemia",
    "lymphoma", "myeloma",
]
SCHEMA = """
CREATE TABLE corpus_stats (avg_doc_length REAL NOT NULL, doc_count INTEGER NOT NULL);
CREATE TABLE idf (token TEXT PRIMARY KEY, idf REAL NOT NULL);
CREATE TABLE documents (
    pmid TEXT PRIMARY KEY, title TEXT, abstract TEXT, journal TEXT, journal_iso TEXT,
    publication_year INTEGER, doc_length INTEGER, mesh_terms TEXT, publication_types TEXT,
    authors TEXT, article_ids TEXT, publication_date TEXT);
CREATE TABLE postings (token TEXT NOT NULL, pmid TEXT NOT NULL, tf INTEGER NOT NULL);
"""


FLAT = False  # set by --flat: source files are <root>/<name>_cancer_inclusive.db (production release layout)


def src_db(root: Path, name: str) -> Path:
    if FLAT:
        return root / f"{name}_cancer_inclusive.db"
    return root / "data" / "pubmed_processed" / f"{name}_cancer_inclusive_all" / "lexical_index.db"


def ro(path: Path) -> sqlite3.Connection:
    return sqlite3.connect("file:" + str(path).replace("\\", "/") + "?mode=ro", uri=True)


def idf_value(n: int, df: int) -> float:
    return math.log(1 + (n - df + 0.5) / (df + 0.5))


def verify_one(name: str, root: Path) -> dict:
    c = ro(src_db(root, name))
    n, avg = c.execute("SELECT doc_count, avg_doc_length FROM corpus_stats").fetchone()
    n_docs = c.execute("SELECT COUNT(*) FROM documents").fetchone()[0]
    mean_len = c.execute("SELECT AVG(doc_length) FROM documents").fetchone()[0]
    null_len = c.execute("SELECT COUNT(*) FROM documents WHERE doc_length IS NULL").fetchone()[0]
    dup = c.execute("SELECT COUNT(*) FROM (SELECT token, pmid FROM postings GROUP BY token, pmid HAVING COUNT(*)>1)").fetchone()[0]
    orphan = c.execute("SELECT COUNT(*) FROM postings p LEFT JOIN documents d ON d.pmid=p.pmid WHERE d.pmid IS NULL").fetchone()[0]
    stored = dict(c.execute("SELECT token, idf FROM idf"))
    maxdiff, mism, missing = 0.0, 0, 0
    seen = set()
    for token, df in c.execute("SELECT token, COUNT(*) FROM postings GROUP BY token"):
        seen.add(token)
        s = stored.get(token)
        if s is None:
            missing += 1
            continue
        d = abs(s - idf_value(n_docs, df))
        maxdiff = max(maxdiff, d)
        if d > 1e-9:
            mism += 1
    extra = len(set(stored) - seen)
    c.close()
    return dict(corpus=name, stored_doc_count=n, actual_docs=n_docs, stored_avg=avg, recomputed_avg=mean_len,
                avg_abs_diff=abs(avg - mean_len), null_doc_length=null_len, duplicate_postings=dup,
                orphan_postings=orphan, idf_tokens_stored=len(stored), idf_tokens_in_postings=len(seen),
                tokens_missing_from_idf=missing, idf_tokens_without_postings=extra,
                idf_max_abs_diff=maxdiff, idf_mismatches=mism)


def build_one(name: str, root: Path, out_root: Path, year_max: int) -> dict:
    t0 = time.time()
    out = out_root / "data" / "pubmed_processed" / f"{name}_cancer_inclusive_all" / "lexical_index.db"
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        out.unlink()
    con = sqlite3.connect("file:" + str(out).replace("\\", "/") + "?mode=rwc", uri=True)
    con.executescript(SCHEMA)
    con.execute("PRAGMA synchronous=OFF")
    con.execute("PRAGMA journal_mode=OFF")
    con.execute("ATTACH DATABASE ? AS src", ("file:" + str(src_db(root, name)).replace("\\", "/") + "?mode=ro",))
    total = con.execute("SELECT COUNT(*) FROM src.documents").fetchone()[0]
    null_year = con.execute("SELECT COUNT(*) FROM src.documents WHERE publication_year IS NULL").fetchone()[0]
    con.execute("INSERT INTO documents SELECT * FROM src.documents WHERE publication_year IS NOT NULL AND publication_year <= ?",
                (year_max,))
    n = con.execute("SELECT COUNT(*) FROM documents").fetchone()[0]
    con.execute("INSERT INTO postings SELECT p.token, p.pmid, p.tf FROM src.postings p "
                "JOIN documents d ON d.pmid = p.pmid")
    con.execute("CREATE INDEX idx_postings_token ON postings(token)")
    avg = con.execute("SELECT AVG(doc_length) FROM documents").fetchone()[0]
    con.execute("INSERT INTO corpus_stats VALUES (?, ?)", (avg, n))
    rows = [(t, idf_value(n, df)) for t, df in con.execute("SELECT token, COUNT(*) FROM postings GROUP BY token")]
    con.executemany("INSERT INTO idf VALUES (?, ?)", rows)
    max_year = con.execute("SELECT MAX(publication_year) FROM documents").fetchone()[0]
    npost = con.execute("SELECT COUNT(*) FROM postings").fetchone()[0]
    con.commit()
    con.execute("DETACH DATABASE src")
    con.close()
    return dict(corpus=name, source_docs=total, source_null_year=null_year, kept_docs=n, dropped_docs=total - n - null_year,
                max_year_kept=max_year, postings=npost, idf_tokens=len(rows), avg_doc_length=avg,
                seconds=round(time.time() - t0))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src-root", type=Path, required=True, help="OncoDigger application root (holds data/pubmed_processed)")
    ap.add_argument("--out-root", type=Path)
    ap.add_argument("--year-max", type=int)
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--flat", action="store_true",
                    help="source layout is <src-root>/<corpus>_cancer_inclusive.db (production release files)")
    ap.add_argument("--corpora", nargs="+", default=CORPORA)
    a = ap.parse_args()
    global FLAT
    FLAT = a.flat
    if a.verify:
        res = [verify_one(n, a.src_root) for n in a.corpora]
        for r in res:
            print(r, flush=True)
        bad = [r for r in res if r["idf_mismatches"] or r["tokens_missing_from_idf"] or r["avg_abs_diff"] > 1e-6
               or r["stored_doc_count"] != r["actual_docs"] or r["orphan_postings"] or r["duplicate_postings"]]
        print("VERIFY", "FAILED for " + ",".join(r["corpus"] for r in bad) if bad else "PASSED for all corpora")
        sys.exit(1 if bad else 0)
    if not (a.out_root and a.year_max):
        ap.error("--out-root and --year-max are required unless --verify")
    res = []
    for n in a.corpora:
        r = build_one(n, a.src_root, a.out_root, a.year_max)
        print(r, flush=True)
        res.append(r)
    rep = a.out_root / f"build_report_year_max_{a.year_max}.csv"
    with rep.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(res[0]))
        w.writeheader()
        w.writerows(res)
    print("wrote", rep)


if __name__ == "__main__":
    main()
