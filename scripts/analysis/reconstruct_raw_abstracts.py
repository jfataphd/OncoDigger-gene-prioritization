"""
reconstruct_raw_abstracts.py -- rebuild and verify the raw abstract files for the 2026-09-27 corpus build.

The audit scripts that read case-preserved title/abstract text (vocab_collision_*, vocab_corrected_rerun*,
met_disambiguation_check, title_case_blind_spot_check, ...) need per-corpus raw JSONL files
(`pubmed_raw/<corpus>/abstracts.jsonl`). The 2026-09-27 build is the 2026-05-08 baseline plus an incremental
PubMed fetch, so no single raw file per corpus exists for it. This script reconstructs one:

  raw(corpus) = records of (May raw UNION incremental raw) whose PMID is in the build's database,
                taking the incremental record where a PMID occurs in both (the later PubMed version).

and then verifies every reconstructed record against the database row: title and abstract are compared
after the same normalisation the index applies (lower-casing, accent folding, HTML-entity decoding, removal of
all non-alphanumerics). Any record that still differs is listed in the manifest.

Environment / arguments (nothing is read from a hard-coded location):
  ONCODIGGER_ROOT   root holding data/pubmed_processed/<corpus>_all/lexical_index.db (the 2026-09-27 build)
  MAY_RAW_DIR       directory holding <corpus>/abstracts.jsonl for the 2026-05-08 baseline
  INCREMENTAL_RAW_DIR   directory holding <corpus>/abstracts.jsonl for the incremental fetch
Output: ONCODIGGER_ROOT/data/pubmed_raw/<corpus>/abstracts.jsonl and data/provenance/raw_reconstruction_2026-09-27.json
"""
from __future__ import annotations

import hashlib
import html
import json
import os
import re
import sqlite3
import unicodedata
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NAMES = ["breast", "lung", "colorectal", "prostate", "melanoma", "bladder", "kidney", "pancreas", "liver", "stomach",
         "esophageal", "ovarian", "endometrial", "cervical", "thyroid", "brain", "leukemia", "lymphoma", "myeloma"]


def key(s: str) -> str:
    s = "".join(c for c in unicodedata.normalize("NFKD", s or "") if not unicodedata.combining(c))
    return re.sub(r"[^a-z0-9]+", "", s.lower())


def load(path: Path) -> dict:
    out = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.strip():
                r = json.loads(line)
                out[str(r["pmid"])] = r
    return out


def main() -> None:
    root = Path(os.environ["ONCODIGGER_ROOT"])
    may, inc = Path(os.environ["MAY_RAW_DIR"]), Path(os.environ["INCREMENTAL_RAW_DIR"])
    manifest, total_bad = {}, 0
    for n in NAMES:
        db = sqlite3.connect(f"file:{(root / 'data' / 'pubmed_processed' / f'{n}_cancer_inclusive_all' / 'lexical_index.db').as_posix()}?mode=ro&immutable=1", uri=True)
        rows = {str(p): (t or "", a or "") for p, t, a in db.execute("select pmid,title,abstract from documents")}
        db.close()
        m, i = load(may / f"{n}_cancer_inclusive" / "abstracts.jsonl"), load(inc / f"{n}_cancer_inclusive" / "abstracts.jsonl")
        missing = set(rows) - (set(m) | set(i))
        assert not missing, (n, len(missing))
        out_dir = root / "data" / "pubmed_raw" / f"{n}_cancer_inclusive"
        out_dir.mkdir(parents=True, exist_ok=True)
        h = hashlib.sha256()
        n_inc = n_may = 0
        mismatched = []
        with open(out_dir / "abstracts.jsonl", "w", encoding="utf-8", newline="\n") as f:
            for pmid in sorted(rows, key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x)):
                rec = i.get(pmid) or m.get(pmid)
                n_inc += pmid in i
                n_may += pmid not in i
                t, a = rows[pmid]
                ok = all(key(rec.get(fld)) == key(db_val) or key(html.unescape(rec.get(fld) or "")) == key(db_val)
                         for fld, db_val in (("title", t), ("abstract", a)))
                if not ok:
                    mismatched.append(pmid)
                line = json.dumps(rec, ensure_ascii=False) + "\n"
                h.update(line.encode("utf-8"))
                f.write(line)
        total_bad += len(mismatched)
        manifest[n] = {"documents": len(rows), "from_incremental": n_inc, "from_may_baseline": n_may,
                       "records_in_sources_not_in_build": len(set(m) | set(i)) - len(rows),
                       "title_or_abstract_differs_from_database": mismatched, "sha256": h.hexdigest()}
        print(f"{n:11s} docs={len(rows):7d} incremental={n_inc:6d} may={n_may:7d} differing={len(mismatched)}", flush=True)
    out = {"build": "2026-09-27", "rule": "(May raw U incremental raw) restricted to the build database PMIDs; incremental preferred on overlap",
           "total_records": sum(v["documents"] for v in manifest.values()), "total_differing": total_bad, "corpora": manifest}
    dst = ROOT / "data" / "provenance" / "raw_reconstruction_2026-09-27.json"
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("Wrote", dst, "| total differing records:", total_bad)


if __name__ == "__main__":
    main()
