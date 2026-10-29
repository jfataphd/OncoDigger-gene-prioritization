"""
evidence_audit_worksheet.py -- build the review worksheet for the evidence behind the September-27 rankings that
were NOT covered by the earlier (May-text) manual checks.

Coverage on the final 2026-09-27 rankings:
  * NEW PAIRS: the four gene-cancer pairs in the final top 10 that were not in the top 10 of the submitted
    (v1.2, May-build) ranking: MTOR (Cervical), PIK3CA (Melanoma), MEN1 (Pancreatic), PALB2 (Prostate).
    20 randomly sampled supporting abstracts each (seed 42), with the sentence around the symbol.
  * SENSITIVE GENES: MET (Lung) and FH (Kidney), after the frozen correction: every retained posting from a
    record with publication year 2026 (the months the May baseline lacks) plus a random 20 of the older ones.
  * RANDOM BLIND SAMPLE: 2 supporting abstracts for each of 30 randomly chosen other top-10 pairs.
The worksheet leaves the `genuine_gene_mention` and `supports_pair` columns empty for the reader.

Historical audit: its pair list is the pre-correction top 10 (data/derived/evidence_audit/audited_top10_pairs.csv).
Env: ONCODIGGER_ROOT = verified 2026-09-27 corpus root. Output: data/derived/evidence_audit/worksheet.csv
Run: python -X utf8 scripts/analysis/evidence_audit_worksheet.py
"""
from __future__ import annotations

import os
import random
import re
import sqlite3
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CORP = Path(os.environ["ONCODIGGER_ROOT"]) / "data" / "pubmed_processed"
CORR = ROOT / "data" / "derived" / "intogen_rescoring" / "corrected_corpora"
OUT = ROOT / "data" / "derived" / "evidence_audit"
KEY = {"Breast Cancer": "breast", "Lung Cancer": "lung", "Colorectal Cancer": "colorectal", "Prostate Cancer": "prostate",
       "Melanoma": "melanoma", "Bladder Cancer": "bladder", "Kidney Cancer": "kidney", "Pancreatic Cancer": "pancreas",
       "Liver Cancer": "liver", "Stomach Cancer": "stomach", "Esophageal Cancer": "esophageal", "Ovarian Cancer": "ovarian",
       "Endometrial Cancer": "endometrial", "Cervical Cancer": "cervical", "Thyroid Cancer": "thyroid", "Brain Cancer": "brain",
       "Leukemia": "leukemia", "Lymphoma": "lymphoma", "Myeloma": "myeloma"}
NEW = [("Cervical Cancer", "MTOR"), ("Melanoma", "PIK3CA"), ("Pancreatic Cancer", "MEN1"), ("Prostate Cancer", "PALB2")]
SENSITIVE = [("Lung Cancer", "MET"), ("Kidney Cancer", "FH")]


def db(cancer: str, corrected: bool) -> sqlite3.Connection:
    k = KEY[cancer]
    p = (CORR if corrected else CORP) / f"{k}_cancer_inclusive_all" / "lexical_index.db"
    return sqlite3.connect(f"file:{p.as_posix()}?mode=ro&immutable=1", uri=True)


def rows(cancer: str, gene: str, kind: str, pmids: list[str], corrected: bool) -> list[dict]:
    con = db(cancer, corrected)
    pat = re.compile(r".{0,120}\b" + re.escape(gene) + r"\b.{0,120}", re.IGNORECASE)
    out = []
    for p in pmids:
        t, a, y = con.execute("SELECT title, abstract, publication_year FROM documents WHERE pmid = ?", (p,)).fetchone()
        txt = f"{t or ''} . {a or ''}".replace("\n", " ")
        m = pat.search(txt)
        out.append(dict(kind=kind, cancer=cancer, gene=gene, pmid=p, year=y, title=t,
                        context=m.group(0) if m else "", genuine_gene_mention="", supports_pair="", note=""))
    con.close()
    return out


def postings(cancer: str, gene: str, corrected: bool) -> list[str]:
    con = db(cancer, corrected)
    r = sorted(str(x[0]) for x in con.execute("SELECT DISTINCT pmid FROM postings WHERE token = ?", (gene.lower(),)))
    con.close()
    return r


def main() -> None:
    rng = random.Random(42)
    # The audit was drawn from the top-10 pairs of the 2026-09-27 ranking as it stood BEFORE the exact-background-count correction
    # (commit 647df5a). Those 190 pairs are kept in audited_top10_pairs.csv so the worksheet stays reproducible; the final
    # ranking in data/canonical differs (for example Melanoma/PIK3CA is no longer in its top 10).
    top = pd.read_csv(OUT / "audited_top10_pairs.csv")
    out: list[dict] = []
    for c, g in NEW:
        assert ((top.cancer == c) & (top.gene == g)).any(), (c, g)
        p = postings(c, g, False)
        out += rows(c, g, "new_pair", rng.sample(p, min(20, len(p))), False)
    for c, g in SENSITIVE:
        assert ((top.cancer == c) & (top.gene == g)).any(), (c, g)
        p = postings(c, g, True)
        con = db(c, True)
        new26 = [x for x in p if (con.execute("SELECT publication_year FROM documents WHERE pmid = ?", (x,)).fetchone() or [0])[0] == 2026]
        con.close()
        old = [x for x in p if x not in set(new26)]
        out += rows(c, g, "sensitive_2026_records", new26, True)
        out += rows(c, g, "sensitive_random_older", rng.sample(old, min(20, len(old))), True)
    done = set(NEW) | set(SENSITIVE)
    others = [(r.cancer, r.gene) for r in top.itertuples() if (r.cancer, r.gene) not in done]
    for c, g in rng.sample(others, 30):
        p = postings(c, g, False)
        out += rows(c, g, "random_blind", rng.sample(p, min(2, len(p))), False)
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(out)
    df.to_csv(OUT / "worksheet.csv", index=False)
    print(df.groupby(["kind", "cancer", "gene"]).size().to_string())
    print("wrote", OUT / "worksheet.csv")


if __name__ == "__main__":
    main()
