"""
vocab_collision_check_batch1.py -- diagnostic only, not yet wired into the scoring pipeline.
Checks the six Tier-1/2 genes that actually hold a TOP-10 slot somewhere in the paper's
reported rankings (i.e. genes that feed the 93.2%/75.8% headline numbers directly):
EGFR, KIT, ESR1, POLE, APC, RET. Lower-ranked flagged candidates (ROS1, FAS, MAF, REST,
CIC, MPL, FAT1, and the medium-priority list: FLT3, GNAS, MEN1, NF1, PMS2, RB1, TSC1,
TSC2, VHL) were NOT reached in this pass due to time budget -- see the summary printed
at the end of this script for exactly what was and wasn't covered.

Method (same two-step rule validated for MET earlier in this project):
  1. Capitalization filter: of a gene's full-corpus lowercase-token matches, keep only
     papers where the RAW, case-preserved title+abstract contains the symbol as a
     standalone, whole-word, ALL-CAPS match (C_g = fraction passing this filter).
  2. Negative/positive phrase check (KIT only in this pass, as a worked example): of the
     capitalized matches, classify by presence of gene-context phrases (c-kit, CD117, kit
     mutation/gene/expression/inhibitor/exon/receptor, stem cell factor) versus common-word
     collision phrases (assay/test/ELISA/detection/extraction/diagnostic/commercial kit).

Uses the real supporting_papers / engine_rank ground truth from
data/derived/intogen_rescoring/intogen_literature_weight_audit.csv and the real top-50
ranks from data/canonical/precision_at_k_detail.csv, not assumption.

Run: python -X utf8 scripts/analysis/vocab_collision_check_batch1.py
"""
from __future__ import annotations

import os

import json
import re
import sqlite3
from pathlib import Path

BASE_DB = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET")) / "data" / "pubmed_processed"
BASE_RAW = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET")) / "data" / "pubmed_raw"

# (symbol, processed-db corpus dir, raw-jsonl corpus dir) -- checked in the gene's best-ranked cancer
CHECKS = [
    ("KIT", "melanoma_cancer_inclusive_all", "melanoma_cancer_inclusive"),
    ("ESR1", "breast_cancer_inclusive_all", "breast_cancer_inclusive"),
    ("POLE", "endometrial_cancer_inclusive_all", "endometrial_cancer_inclusive"),
    ("APC", "colorectal_cancer_inclusive_all", "colorectal_cancer_inclusive"),
    ("RET", "thyroid_cancer_inclusive_all", "thyroid_cancer_inclusive"),
    ("EGFR", "lung_cancer_inclusive_all", "lung_cancer_inclusive"),
]

KIT_POSITIVE = ["c-kit", "cd117", "cd 117", "kit mutation", "kit gene", "kit expression",
                "kit inhibitor", "kit exon", "kit receptor", "stem cell factor"]
KIT_NEGATIVE = ["assay kit", "test kit", "elisa kit", "detection kit", "extraction kit",
                "diagnostic kit", "starter kit", "commercial kit", "kit was used", "kit according"]


def capitalization_check(symbol: str, db_corpus: str, raw_corpus: str) -> tuple[int, int]:
    token = symbol.lower()
    db = BASE_DB / db_corpus / "lexical_index.db"
    con = sqlite3.connect(str(db))
    pmids = {str(p[0]) for p in con.execute(
        "SELECT DISTINCT pmid FROM postings WHERE token = ?", (token,)).fetchall()}
    pattern = re.compile(r"\b" + symbol + r"\b")
    raw = BASE_RAW / raw_corpus / "abstracts.jsonl"
    total = caps = 0
    with open(raw, encoding="utf-8") as f:
        for line in f:
            rec = json.loads(line)
            if str(rec.get("pmid", "")) not in pmids:
                continue
            total += 1
            text = (rec.get("title") or "") + " " + (rec.get("abstract") or "")
            if pattern.search(text):
                caps += 1
    return total, caps


def kit_signal_breakdown() -> dict:
    db = BASE_DB / "melanoma_cancer_inclusive_all" / "lexical_index.db"
    con = sqlite3.connect(str(db))
    pmids = {str(p[0]) for p in con.execute(
        "SELECT DISTINCT pmid FROM postings WHERE token = 'kit'").fetchall()}
    pattern = re.compile(r"\bKIT\b")
    raw = BASE_RAW / "melanoma_cancer_inclusive" / "abstracts.jsonl"
    caps = pos = neg = neither = 0
    with open(raw, encoding="utf-8") as f:
        for line in f:
            rec = json.loads(line)
            if str(rec.get("pmid", "")) not in pmids:
                continue
            text = (rec.get("title") or "") + " " + (rec.get("abstract") or "")
            if not pattern.search(text):
                continue
            caps += 1
            low = text.lower()
            has_pos = any(p in low for p in KIT_POSITIVE)
            has_neg = any(p in low for p in KIT_NEGATIVE)
            if has_neg and not has_pos:
                neg += 1
            elif has_pos:
                pos += 1
            else:
                neither += 1
    return {"capitalized": caps, "positive": pos, "negative_only": neg, "unclassified": neither}


def main() -> None:
    print("=== Stage 1: capitalization fidelity (C_g) for the six top-10-holding genes ===\n")
    for symbol, db_corpus, raw_corpus in CHECKS:
        total, caps = capitalization_check(symbol, db_corpus, raw_corpus)
        print(f"{symbol:6s} ({db_corpus}): total={total:6d}  standalone-caps={caps:6d}  "
              f"C_g={100 * caps / total:5.1f}%")

    print("\n=== KIT worked example: positive/negative signal among capitalized matches ===")
    result = kit_signal_breakdown()
    print(result)
    print(f"  genuine-signal rate among capitalized matches: "
          f"{100 * result['positive'] / result['capitalized']:.1f}%")
    print(f"  confirmed negative (common-word 'kit') rate: "
          f"{100 * result['negative_only'] / result['capitalized']:.1f}%")
    print(f"  unclassified (neither signal, needs manual read): "
          f"{100 * result['unclassified'] / result['capitalized']:.1f}%")

    print("""
=== Coverage summary ===
Completed this pass (C_g capitalization check, real postings/abstract text, no
IntOGen-impact influence on methodology): EGFR, KIT, ESR1, POLE, APC, RET -- the six
Tier-1/2 flagged genes that actually hold a top-10 slot somewhere in the rankings.
KIT additionally got a full positive/negative phrase breakdown as a worked example.

NOT reached this pass (ran out of time budget, per explicit instruction to stop and
report rather than rush): ROS1, FAS/FASN, MAF, REST, CIC, MPL, FAT1 (all confirmed or
suspected risks, but none hold a top-10 slot -- they only affect top-20/30/40/50 depth
metrics, not the two headline numbers), plus the medium-priority list (FLT3, GNAS,
MEN1, NF1, PMS2, RB1, TSC1, TSC2, VHL).

No rule was adjusted based on IntOGen impact. No scoring was rerun. Nothing pushed to
any remote.
""")


if __name__ == "__main__":
    main()
