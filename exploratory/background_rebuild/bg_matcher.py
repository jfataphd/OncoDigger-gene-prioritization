"""
bg_matcher.py -- the engine's gene matcher applied to whole collections, with version-stamped, atomically written counts.
(Implements docs/BACKGROUND_COUNT_SPEC.md. Exploratory until the correction is frozen.)

  build(root)                    -> Matcher with the SAME lookups the paper pipeline passes to associated_entity_tables
  Matcher.entities_in(doc)       -> {(entity_id, label)} found in one paper (copy of the matching block of associated_entity_tables)
  Matcher.count(rows, pool_ids)  -> (B per entity, S per entity within pool_ids, papers); each PMID counted once
  stamp(db_path, root)           -> dict that identifies corpus content, vocabulary, matching policy, redaction policy and matcher version
  save(path_prefix, ...) / load(path_prefix, expected_stamp) -> atomic write; load refuses any stamp mismatch
"""
import hashlib
import importlib.util
import inspect
import json
import os
import sys
import tempfile
from collections import Counter
from pathlib import Path

MATCHER_VERSION = "bgmatch-1"
REDACTION_POLICY = "MET, REST and FH false-positive mentions redacted from title and abstract text of the corrected databases (Methods 3.3)"
REPO = Path(__file__).resolve().parents[2]


def _sha256(path, block=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


class Matcher:
    def __init__(self, root: Path, query_text: str | None = None):
        sys.path.insert(0, str(REPO / "scripts" / "analysis"))
        spec = importlib.util.spec_from_file_location("rqp", REPO / "scripts" / "analysis" / "run_query_planes.py")
        self.m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.m)
        sys.path.insert(0, str(root / "src"))
        from oncodigger.search import discovery
        from oncodigger.search.lexical import normalize_query
        from oncodigger.search.lexicons import load_lexicons
        self.root, self.discovery = Path(root), discovery
        lex = load_lexicons(self.root)
        self.lookup = self.m.longest_match_lookups(self.m.suppress_ambiguous_symbols(lex["Genes"]), discovery)
        tokens = normalize_query(query_text if query_text is not None else self.m.PLANES["mutation"])
        self.query_terms = {t.replace("-", " ") for t in tokens}
        self.single, self.phrases = self.lookup["single"], self.lookup["phrases"]

    def entities_in(self, doc: dict) -> set:
        d = self.discovery
        text = d.full_document_text(doc)
        padded = f" {text} "
        token_set = set(d.document_tokens(doc))
        seen = set()
        for token in token_set:
            if token in self.query_terms:
                continue
            for entity in self.single.get(token, ()):
                seen.add((entity["id"], entity["label"]))
        for anchor in token_set:
            for phrase, entity in self.phrases.get(anchor, ()):
                if phrase in self.query_terms:
                    continue
                if f" {phrase} " in padded:
                    seen.add((entity["id"], entity["label"]))
        return {k for k in seen if d.norm_text(k[1]) not in self.query_terms}

    def count(self, rows, pool_ids=frozenset(), progress=None):
        """rows: iterable of (pmid, title, abstract). A PMID that occurs twice is counted once."""
        b, s, done, n = Counter(), Counter(), set(), 0
        for pmid, title, abstract in rows:
            pmid = str(pmid)
            if pmid in done:
                continue
            done.add(pmid)
            found = self.entities_in({"title": title or "", "abstract": abstract or ""})
            inpool = pmid in pool_ids
            for k in found:
                b[k] += 1
                if inpool:
                    s[k] += 1
            n += 1
            if progress and n % 25000 == 0:
                progress(n)
        return b, s, n

    def policy_hash(self) -> str:
        src = "".join(inspect.getsource(f) for f in (self.m.longest_match_lookups, self.m.suppress_ambiguous_symbols))
        src += repr(sorted(self.m.SUPPRESS_SYMBOLS)) + repr(self.m.GENE_REMAP) + repr(sorted(self.m.GENE_EXCLUDE))
        return hashlib.sha256(src.encode()).hexdigest()


def stamp(matcher: Matcher, db_path: Path) -> dict:
    root = matcher.root
    code = {}
    for rel in ("src/oncodigger/search/discovery.py", "src/oncodigger/search/lexicons.py", "src/oncodigger/search/association.py"):
        code[rel] = _sha256(root / rel)
    return {"matcher_version": MATCHER_VERSION, "corpus_db": Path(db_path).parent.name, "corpus_db_sha256": _sha256(db_path),
            "vocabulary_csv_sha256": _sha256(root / "Human Genes.csv"), "matching_policy_sha256": matcher.policy_hash(),
            "engine_code_sha256": code, "redaction_policy": REDACTION_POLICY}


def save(prefix: Path, header: list, rows: list, st: dict):
    """Atomic: write both files to temporary names in the same directory, then rename (stamp last, so a stamp implies a complete table)."""
    prefix = Path(prefix)
    d = prefix.parent
    fd, tmp_csv = tempfile.mkstemp(dir=d, suffix=".tmp")
    with os.fdopen(fd, "w", newline="", encoding="utf-8") as fh:
        import csv
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    os.replace(tmp_csv, f"{prefix}.csv")
    fd, tmp_js = tempfile.mkstemp(dir=d, suffix=".tmp")
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        json.dump(st, fh, indent=1, sort_keys=True)
    os.replace(tmp_js, f"{prefix}.stamp.json")


def load(prefix: Path, expected: dict):
    import pandas as pd
    prefix = Path(prefix)
    got = json.load(open(f"{prefix}.stamp.json", encoding="utf-8"))
    if got != expected:
        diff = sorted(k for k in set(got) | set(expected) if got.get(k) != expected.get(k))
        raise RuntimeError(f"stale background counts for {prefix.name}: stamp differs in {diff}")
    return pd.read_csv(f"{prefix}.csv")
