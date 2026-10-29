"""
export_pmid_manifests.py -- the PMID list of each of the 19 collections of the 2026-09-27 build, with counts and SHA-256, so the exact set of
records behind every result is identified without redistributing the databases (about 27 GB).

Output: data/provenance/pmid_manifests/<collection>.pmids.gz (one PMID per line, ascending; gzip with mtime 0 so the file is reproducible),
        data/provenance/pmid_manifests/MANIFEST.csv (collection, n_pmids, sha256 of the decompressed list, database sha256 from
        corpus_build_2026-09-27.json), and the pooled (de-duplicated) count.
Run:    ONCODIGGER_ROOT=<root with data/pubmed_processed> python -X utf8 scripts/reproduce/export_pmid_manifests.py [--check]
--check recomputes the lists from the databases and verifies the committed files without writing.
"""
import csv
import gzip
import hashlib
import json
import os
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "provenance" / "pmid_manifests"
NAMES = ["breast", "lung", "colorectal", "prostate", "melanoma", "bladder", "kidney", "pancreas", "liver", "stomach", "esophageal",
         "ovarian", "endometrial", "cervical", "thyroid", "brain", "leukemia", "lymphoma", "myeloma"]


def main() -> None:
    check = "--check" in sys.argv
    base = Path(os.environ["ONCODIGGER_ROOT"]) / "data" / "pubmed_processed"
    build = json.loads((ROOT / "data" / "provenance" / "corpus_build_2026-09-27.json").read_text(encoding="utf-8"))
    dbsha = {d["corpus"]: d["sha256"] for d in build["corpora"]}
    OUT.mkdir(parents=True, exist_ok=True)
    rows, allp, bad = [], set(), 0
    for k in NAMES:
        con = sqlite3.connect(f"file:{(base / f'{k}_cancer_inclusive_all' / 'lexical_index.db').as_posix()}?mode=ro&immutable=1", uri=True)
        pm = sorted((str(r[0]) for r in con.execute("select pmid from documents")), key=int)
        con.close()
        text = ("\n".join(pm) + "\n").encode()
        sha = hashlib.sha256(text).hexdigest()
        allp.update(pm)
        f = OUT / f"{k}.pmids.gz"
        if check:
            have = gzip.decompress(f.read_bytes())
            ok = have == text
            bad += not ok
            print(f"{k}: {'OK' if ok else 'MISMATCH'} ({len(pm)} PMIDs)")
        else:
            with f.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as g:
                g.write(text)
        rows.append({"collection": k, "n_pmids": len(pm), "pmid_list_sha256": sha, "database_sha256": dbsha.get(k, "")})
    if check:
        sys.exit(1 if bad else 0)
    with (OUT / "MANIFEST.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} lists; {sum(r['n_pmids'] for r in rows)} collection records, {len(allp)} distinct PMIDs")


if __name__ == "__main__":
    main()
