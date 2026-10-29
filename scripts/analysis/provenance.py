"""
provenance.py -- guard that a corpus root really is the recorded corpus build.

Why this exists: earlier outputs of this project were produced from two different corpus builds (the
2026-05-08 baseline and the 2026-09-27 refresh) without anything stopping or flagging it, because the engine
scripts silently fell back to a default corpus folder. This module makes the build explicit.

assert_corpus_build(root) reads data/provenance/corpus_build_<BUILD_ID>.json and checks, for all 19 corpora
under `root/data/pubmed_processed/<corpus>_all/lexical_index.db`:
  * the database exists, its file size equals the recorded size, and
  * its document count (documents table and corpus_stats) equals the recorded count.
With full=True it also recomputes each file's SHA-256 (slow: the 19 databases total about 27 GB).
It raises RuntimeError (never silently continues) when anything differs. Scripts that read the corpus
databases call it at start-up (run_query_planes.py does).
"""
from __future__ import annotations

import hashlib
import json
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILD_ID = "2026-09-27"
MANIFEST = ROOT / "data" / "provenance" / f"corpus_build_{BUILD_ID}.json"


def load_manifest() -> dict:
    return json.loads(MANIFEST.read_text(encoding="utf-8"))


TEXT_SUFFIXES = {".csv", ".json", ".md", ".py", ".txt", ".bak", ".tsv"}


def file_digest(path) -> str:
    """SHA-256 of a committed output. Text files are hashed with CRLF normalised to LF so the digest is the same in a
    Windows checkout (core.autocrlf) and on Linux; binary files are hashed byte for byte."""
    path = Path(path)
    data = path.read_bytes()
    if path.suffix.lower() in TEXT_SUFFIXES:
        data = data.replace(b"\r\n", b"\n")
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def assert_corpus_build(root, full: bool = False, verbose: bool = True) -> dict:
    root = Path(root)
    man = load_manifest()
    problems = []
    for c in man["corpora"]:
        db = root / "data" / "pubmed_processed" / f"{c['corpus']}_cancer_inclusive_all" / "lexical_index.db"
        if not db.exists():
            problems.append(f"{c['corpus']}: database not found at {db}")
            continue
        if db.stat().st_size != c["bytes"]:
            problems.append(f"{c['corpus']}: file size {db.stat().st_size} != recorded {c['bytes']}")
        con = sqlite3.connect(f"file:{db.as_posix()}?mode=ro&immutable=1", uri=True)
        n_docs = con.execute("select count(*) from documents").fetchone()[0]
        stored = con.execute("select doc_count from corpus_stats").fetchone()[0]
        con.close()
        if n_docs != c["documents"] or stored != c["documents"]:
            problems.append(f"{c['corpus']}: documents {n_docs} (stored {stored}) != recorded {c['documents']}")
        if full and sha256_file(db) != c["sha256"]:
            problems.append(f"{c['corpus']}: SHA-256 differs from the recorded build")
    if problems:
        raise RuntimeError(
            f"Corpus root {root} is NOT the recorded {BUILD_ID} build:\n  " + "\n  ".join(problems)
            + "\nSet ONCODIGGER_ROOT to the verified 2026-09-27 corpus folder (see data/provenance/README.md).")
    if verbose:
        total = sum(c["documents"] for c in man["corpora"])
        print(f"[provenance] corpus root verified as build {BUILD_ID}: 19 databases, {total:,} documents"
              f"{' (SHA-256 checked)' if full else ' (size + document count checked)'}", flush=True)
    return man
