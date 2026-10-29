"""
cutoff_matched_vocab_correction.py -- apply the MET/REST/FH correction to the date-truncated corpora.

The main analysis reports the ranking after the frozen MET (Lung Cancer), REST (Thyroid Cancer) and FH
(Kidney Cancer) false-positive redaction (Methods 3.3). The cutoff-matched sensitivity analysis (Methods 2.7)
must use the same policy so that publication year is the only variable that changes. This script:

  1. identifies the false-positive PMIDs of each gene on the full 2026-09-27 build with the frozen rule
     (vocab_corrected_rerun.false_positive_pmids, unchanged); the rule is per posting, so its decision for a PMID
     does not depend on the corpus around it;
  2. for each cutoff, copies the truncated Lung / Thyroid / Kidney database, redacts the gene token from the
     title and abstract of every false-positive PMID that the truncated database contains, and deletes the
     gene's postings rows for those PMIDs (as make_corrected_db does on the full build; IDF and corpus
     statistics are left as in the truncated database);
  3. re-ranks those three cancers with the real engine (run_query_planes.rank_genes, mutation plane);
  4. writes data/derived/cutoff_matched_rankings/mutation_corpus_to_<year>_vocab_corrected.csv = the uncorrected
     truncated ranking with those three cancers replaced, plus a report of how many false-positive PMIDs fell
     inside each truncated corpus.

Env: ONCODIGGER_ROOT = verified 2026-09-27 corpus root (engine source, lexicons, raw abstracts).
Args: --trunc-root-2022 / --trunc-root-2023 = output roots of build_date_truncated_corpora.py;
      --work-dir = scratch folder for the corrected copies (several GB; not kept in the repository).
Run: PYTHONHASHSEED=0 python -X utf8 scripts/analysis/cutoff_matched_vocab_correction.py ...
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

TARGETS = [("MET", "Lung Cancer", "lung_cancer_inclusive"),
           ("REST", "Thyroid Cancer", "thyroid_cancer_inclusive"),
           ("FH", "Kidney Cancer", "kidney_cancer_inclusive")]
PLACEHOLDER = "ZZREDACTEDZZ"
RANK_DIR = ROOT / "data" / "derived" / "cutoff_matched_rankings"


def redact_copy(src: Path, dst: Path, gene: str, fp: set[str]) -> dict:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    con = sqlite3.connect(str(dst))
    token = gene.lower()
    pat = re.compile(chr(92) + "b" + re.escape(token) + chr(92) + "b", re.IGNORECASE)
    present = [p for p in fp if con.execute("SELECT 1 FROM documents WHERE pmid = ?", (p,)).fetchone()]
    for pmid in present:
        title, abstract = con.execute("SELECT title, abstract FROM documents WHERE pmid = ?", (pmid,)).fetchone()
        con.execute("UPDATE documents SET title = ?, abstract = ? WHERE pmid = ?",
                    (pat.sub(PLACEHOLDER, title or ""), pat.sub(PLACEHOLDER, abstract or ""), pmid))
    before = con.execute("SELECT COUNT(*) FROM postings WHERE token = ?", (token,)).fetchone()[0]
    if present:
        con.execute("DELETE FROM postings WHERE token = ? AND pmid IN (" + ",".join("?" * len(present)) + ")",
                    [token, *present])
    con.commit()
    after = con.execute("SELECT COUNT(*) FROM postings WHERE token = ?", (token,)).fetchone()[0]
    con.close()
    return dict(fp_on_full_build=len(fp), fp_in_truncated_corpus=len(present), postings_before=before, postings_after=after)


def main() -> None:
    if os.environ.get("PYTHONHASHSEED") != "0":
        raise SystemExit("Run with PYTHONHASHSEED=0")
    ap = argparse.ArgumentParser()
    ap.add_argument("--trunc-root-2022", type=Path, required=True)
    ap.add_argument("--trunc-root-2023", type=Path, required=True)
    ap.add_argument("--work-dir", type=Path, required=True)
    a = ap.parse_args()
    root = Path(os.environ["ONCODIGGER_ROOT"])

    import vocab_corrected_rerun as vcr
    from run_query_planes import suppress_ambiguous_symbols, longest_match_lookups, open_corpus, rank_genes, PLANES
    sys.path.insert(0, str(root / "src"))
    from oncodigger.search import discovery
    from oncodigger.search.discovery import associated_entity_tables
    from oncodigger.search.lexical import bm25_search, normalize_query
    from oncodigger.search.lexicons import load_lexicons

    fps = {g: vcr.false_positive_pmids(g, c) for g, c, _ in TARGETS}
    lex = load_lexicons(root)
    genes = longest_match_lookups(suppress_ambiguous_symbols(lex["Genes"]), discovery)
    engine = (bm25_search, normalize_query, associated_entity_tables, {"Genes": genes})

    import pandas as pd
    report = []
    for year, troot in ((2022, a.trunc_root_2022), (2023, a.trunc_root_2023)):
        base = pd.read_csv(RANK_DIR / f"mutation_corpus_to_{year}.csv")
        for gene, cancer, key in TARGETS:
            src = troot / "data" / "pubmed_processed" / f"{key}_all" / "lexical_index.db"
            dst = a.work_dir / f"y{year}" / f"{key}_all" / "lexical_index.db"
            # control: the engine on the UNCORRECTED truncated database must reproduce the committed ranking
            con = open_corpus(src, discovery)
            ctrl = rank_genes(con, PLANES["mutation"], engine)
            con.close()
            committed = list(base[base.cancer == cancer].sort_values("rank").gene)
            assert ctrl[:len(committed)] == committed, f"control run differs from committed ranking: {year} {cancer}"
            stats = redact_copy(src, dst, gene, fps[gene])
            con = open_corpus(dst, discovery)
            ranked = rank_genes(con, PLANES["mutation"], engine)
            con.close()
            old = base[base.cancer == cancer].sort_values("rank")
            old_top = list(old.head(50).gene)
            new_rows = pd.DataFrame({"cancer": cancer, "rank": range(1, len(ranked) + 1), "gene": ranked})
            keep = [c for c in base.columns if c in new_rows.columns]
            base = pd.concat([base[base.cancer != cancer], new_rows[keep]], ignore_index=True)
            report.append(dict(cutoff=year, cancer=cancer, gene=gene, **stats,
                               top10_changed=int(list(ranked[:10]) != old_top[:10]),
                               rank_before=(old_top.index(gene) + 1) if gene in old_top else "",
                               rank_after=(ranked.index(gene) + 1) if gene in ranked else ""))
            print(report[-1], flush=True)
        order = {c: i for i, c in enumerate(pd.read_csv(RANK_DIR / f"mutation_corpus_to_{year}.csv").cancer.unique())}
        base["_o"] = base.cancer.map(order)
        base = base.sort_values(["_o", "rank"]).drop(columns="_o")
        base.to_csv(RANK_DIR / f"mutation_corpus_to_{year}_vocab_corrected.csv", index=False)
    with (RANK_DIR / "vocab_correction_report.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(report[0]))
        w.writeheader()
        w.writerows(report)
    print("wrote", RANK_DIR)


if __name__ == "__main__":
    main()
