"""
exact_background_counts.py -- exact background count B for scripts that run their own gene matcher (docs/BACKGROUND_COUNT_SPEC.md).

B(gene) = number of distinct papers of the collection in which the SAME per-paper matching loop that produces the supporting-paper
count S finds the gene. The loop below is the one in rank_vs_enrichment_ablation.compute_raw_terms (single tokens, then phrases by
anchor token and substring test; the query terms and the entity label are excluded), applied to every paper instead of the retrieved
pool. The caller checks that the pool support recomputed here equals its own S for every gene.

exact_background(db_path, vocab_lexicon, query_tokens, discovery_module_name, workers) -> (B per gene id, number of papers)
Papers are split over separate worker processes (rowid % workers); a PMID occurs once in a collection (PRIMARY KEY).
"""
import importlib
import json
import os
import pickle
import sqlite3
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


def match_papers(rows, combined_single, combined_phrases, query_terms, discovery):
    counts = Counter()
    for doc in rows:
        text = discovery.full_document_text(doc)
        padded = f" {text} "
        token_set = set(discovery.document_tokens(doc))
        seen = set()
        for token in token_set:
            if token in query_terms:
                continue
            for entity in combined_single.get(token, []):
                seen.add(entity["id"])
        for anchor in token_set:
            for phrase, entity in combined_phrases.get(anchor, []):
                if phrase in query_terms:
                    continue
                if f" {phrase} " in padded:
                    seen.add(entity["id"])
        for gene in seen:
            counts[gene] += 1
    return counts


def exact_background(db_path, vocab_lexicon, query_tokens, discovery_name, workers=6):
    work = Path(tempfile.mkdtemp(prefix="exactb_"))
    with open(work / "lex.pkl", "wb") as fh:
        pickle.dump(vocab_lexicon, fh, protocol=4)
    procs = []
    for i in range(workers):
        spec = {"sys_path": [p for p in sys.path if p], "db": str(db_path), "lex": str(work / "lex.pkl"), "tokens": list(query_tokens),
                "discovery": discovery_name, "chunk": i, "chunks": workers, "out": str(work / f"p{i}.pkl")}
        (work / f"s{i}.json").write_text(json.dumps(spec))
        env = dict(os.environ, PYTHONHASHSEED="0", PYTHONPATH="")
        env.pop("BGFIX_ENABLE", None)
        procs.append(subprocess.Popen([sys.executable, "-X", "utf8", __file__, str(work / f"s{i}.json")], env=env))
    if any(p.wait() for p in procs):
        raise RuntimeError("exact_background worker failed")
    total, n = Counter(), 0
    for i in range(workers):
        c, k = pickle.load(open(work / f"p{i}.pkl", "rb"))
        total.update(c)
        n += k
    for f in work.iterdir():
        f.unlink()
    work.rmdir()
    return dict(total), n


if __name__ == "__main__":
    spec = json.load(open(sys.argv[1]))
    sys.path[:0] = spec["sys_path"]
    discovery = importlib.import_module(spec["discovery"])
    lex = pickle.load(open(spec["lex"], "rb"))
    single, phrases = defaultdict(list), defaultdict(list)
    for t, es in lex.get("single", {}).items():
        single[t].extend(es)
    for a, pes in lex.get("phrases", {}).items():
        phrases[a].extend(pes)
    qterms = {t.replace("-", " ") for t in spec["tokens"]}
    con = sqlite3.connect("file:" + Path(spec["db"]).as_posix() + "?mode=ro", uri=True)
    cnt = [0]

    def rows():
        for pmid, title, abstract in con.execute("SELECT pmid, title, abstract FROM documents WHERE rowid % ? = ?", (spec["chunks"], spec["chunk"])):
            cnt[0] += 1
            yield {"title": title or "", "abstract": abstract or ""}

    counts = match_papers(rows(), single, phrases, qterms, discovery)
    pickle.dump((counts, cnt[0]), open(spec["out"], "wb"), protocol=4)
