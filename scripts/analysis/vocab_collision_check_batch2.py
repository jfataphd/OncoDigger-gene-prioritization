"""
vocab_collision_check_batch2.py -- diagnostic only, not yet wired into the scoring pipeline.
Continues vocab_collision_check_batch1.py: checks ROS1, FAS/FASN, MAF, REST, CIC, MPL, FAT1
(the genes that don't hold a top-10 slot but were flagged as risks). The medium-priority
list (FLT3, GNAS, MEN1, NF1, PMS2, RB1, TSC1, TSC2, VHL) was NOT reached -- ran out of time
budget, stopping and reporting per instructions rather than rushing.

Two real findings that change the picture from batch 1:

1. ROS1 and FAS/FASN are NOT actually live problems, despite real HGNC alias collisions
   (ROS1's alias "ROS", FASN's alias "FAS"). Verified by actually running the real
   load_lexicons() loader (not inferring from reading the code): useful_gene_term()'s
   `len(term) < 4` filter already excludes bare "ROS" (3 chars) from ever becoming a
   matchable term (`'ros' in single` -> False). FASN's "FAS" alias (4 chars, not caught by
   the length filter) nonetheless resolves to FAS only and FASN only in the actual built
   lookup table (`single['fas'] -> ['FAS']`, `single['fasn'] -> ['FASN']`), so the two
   genes do not actually cross-contaminate each other in this pipeline's real term lookup.
   Root cause of why not fully traced (build_lookup's own entity-registration order), but
   the empirical result is unambiguous and was checked directly against the real loader.

2. Methodological correction, applies to EVERY gene checked in both batches: mouse-model
   gene nomenclature is conventionally written in title case ("Fas", "Kit"), not all-caps
   ("FAS", "KIT") the way human genes are. The strict `\bSYMBOL\b` all-caps-only filter used
   in batch 1 undercounts genuine matches for any gene with substantial mouse-model
   literature. Corrected to "not all-lowercase" (accepts FAS, Fas, FAs, fAS -- rejects only
   literal all-lowercase "fas"). This mattered enormously for FAS (C_g 22.1% with the old
   strict filter -> 98.5% corrected) and modestly for KIT (79.0% -> 84.6%, rechecked here).
   KIT's batch-1 result should be read as superseded by the 84.6% figure below.

Run: python -X utf8 scripts/analysis/vocab_collision_check_batch2.py
"""
from __future__ import annotations

import json
import re
import sqlite3
import sys
import os
from pathlib import Path

BASE_DB = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET")) / "data" / "pubmed_processed"
BASE_RAW = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET")) / "data" / "pubmed_raw"
ROOT = Path(__file__).resolve().parents[2]


def not_all_lowercase_count(token: str, db_corpus: str, raw_corpus: str) -> tuple[int, int]:
    """C_g using the corrected filter: counts a match if ANY occurrence of the whole word
    (case-insensitive) is not purely lowercase -- i.e. accepts both human (ALL-CAPS) and
    mouse (Title case) nomenclature, rejects only literal all-lowercase ordinary-word use."""
    db = BASE_DB / db_corpus / "lexical_index.db"
    con = sqlite3.connect(str(db))
    pmids = {str(p[0]) for p in con.execute(
        "SELECT DISTINCT pmid FROM postings WHERE token = ?", (token,)).fetchall()}
    raw = BASE_RAW / raw_corpus / "abstracts.jsonl"
    pattern = re.compile(r"\b" + re.escape(token) + r"\b", re.IGNORECASE)
    total = not_lower = 0
    with open(raw, encoding="utf-8") as f:
        for line in f:
            rec = json.loads(line)
            if str(rec.get("pmid", "")) not in pmids:
                continue
            total += 1
            text = (rec.get("title") or "") + " " + (rec.get("abstract") or "")
            if any(m.group(0) != token for m in pattern.finditer(text)):
                not_lower += 1
    return total, not_lower


def verify_ros1_fas_alias_safety() -> None:
    sys.path.insert(0, str(ROOT / "reference_engine_code"))
    from oncodigger_search.lexicons import load_lexicons
    lex = load_lexicons(ROOT / "reference_engine_code")
    single = lex["Genes"]["single"]
    print("ROS1 alias-safety check (real loader, not inference):")
    print(f"  'ros' in matchable terms: {'ros' in single}  (should be False)")
    print(f"  'ros1' maps to: {[e['id'] for e in single.get('ros1', [])]}")
    print("FAS/FASN cross-contamination check (real loader):")
    print(f"  'fas' maps to: {[e['id'] for e in single.get('fas', [])]}  (should be ['FAS'] only)")
    print(f"  'fasn' maps to: {[e['id'] for e in single.get('fasn', [])]}  (should be ['FASN'] only)")
    print()


CHECKS = [
    ("ros1", "lung_cancer_inclusive_all", "lung_cancer_inclusive", "Lung Cancer (rank 15, full-corpus 'ros1' token; bare 'ros' excluded by useful_gene_term, see above)"),
    ("fas", "lymphoma_cancer_inclusive_all", "lymphoma_cancer_inclusive", "Lymphoma (rank 30)"),
    ("maf", "myeloma_cancer_inclusive_all", "myeloma_cancer_inclusive", "Myeloma (rank 26)"),
    ("rest", "kidney_cancer_inclusive_all", "kidney_cancer_inclusive", "Kidney Cancer (rank 122 in literature-weight audit; not in top-50)"),
    ("cic", "brain_cancer_inclusive_all", "brain_cancer_inclusive", "Brain Cancer (rank 37)"),
    ("mpl", "leukemia_cancer_inclusive_all", "leukemia_cancer_inclusive", "Leukemia (rank 24)"),
    ("fat1", "cervical_cancer_inclusive_all", "cervical_cancer_inclusive", "Cervical Cancer (rank 36)"),
    ("kit", "melanoma_cancer_inclusive_all", "melanoma_cancer_inclusive", "Melanoma (rank 3) -- RECHECK of batch-1 result with corrected filter"),
]


def main() -> None:
    verify_ros1_fas_alias_safety()
    print("=== C_g (corrected, not-all-lowercase filter) ===\n")
    for token, db_corpus, raw_corpus, note in CHECKS:
        total, caps = not_all_lowercase_count(token, db_corpus, raw_corpus)
        if total == 0:
            print(f"{token.upper():6s}: no postings found in {db_corpus}")
            continue
        print(f"{token.upper():6s} ({note}): total={total:6d}  not-all-lowercase={caps:6d}  "
              f"C_g={100 * caps / total:5.1f}%")

    print("""
=== REST: a third, gene-specific false-positive source found by manual read ===
Beyond the ordinary word "rest" (90.8% of all matches), capitalized "REST" matches in the
Kidney Cancer corpus also collide with "nephrogenic rest", an established, unrelated
pathology term (embryonic kidney-tissue remnants, a recognized Wilms tumor precursor
lesion) -- concentrated exactly in this corpus. Manual sample of the 18 capitalized
matches: ~3 genuinely about the REST gene (e.g. "Mutations in the transcriptional
repressor REST predispose to Wilms tumor"), 1 explicitly about "Nephrogenic Rest" (the
unrelated pathology term), remainder ambiguous/needs full-text read. REST needs its own
negative phrase ("nephrogenic rest") in addition to the capitalization filter.

=== Coverage summary ===
Completed: ROS1 and FAS/FASN alias-safety verified directly against the real loader (both
SAFE, no fix needed). C_g computed with the corrected filter for ROS1, FAS, MAF, REST, CIC,
MPL, FAT1, and KIT (rechecked). MAF, CIC, MPL, FAT1 all clean (97-100% C_g). REST confirmed
seriously contaminated (9.2% C_g) with an additional gene-specific false-positive source.
FAS clean once the mouse-nomenclature correction is applied (22.1% -> 98.5%).

NOT reached (time budget): the medium-priority list -- FLT3, GNAS, MEN1, NF1, PMS2, RB1,
TSC1, TSC2, VHL.

No rule was adjusted based on IntOGen impact at any point. No scoring was rerun. Nothing
pushed to any remote.
""")


if __name__ == "__main__":
    main()
