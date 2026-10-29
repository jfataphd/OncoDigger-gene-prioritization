"""
vocab_inventory_screen.py -- vocabulary-wide screen for alias and symbol collisions (policy v5; Supplementary Methods S1).

Why: the earlier collision screen tested each gene's SYMBOL only and was run on the genes that reached a top-50 list, so it could not see
an alias that matches an ordinary word (RANK, nude, Lobe, race, ...). This screen covers the WHOLE frozen gene vocabulary, whether or not a
gene ever reaches a ranking, using only the collections' token counts (no matcher, no ranking, no reference standard).

Alias arm: every active single-word alias token of every gene (token other than the gene's own symbol). A token is a CANDIDATE when, in at least
one of the 19 collections, it occurs in >= ALIAS_MIN papers and in >= RATIO times as many papers as the gene's own symbol token.
Symbol arm: every active single-word symbol token of 3 or more letters that is an ordinary English word (in the 20,000 most common words,
data/reference/english_common_words_20k.txt, from github.com/first20hours/google-10000-english; not part of this release, the script
fetches it from that repository when it is absent and checks it against the SHA-256 below) and occurs in >= SYMBOL_MIN papers of some collection.
Each candidate is annotated with: the collection with the largest count, the alias and symbol counts there, whether the token is an ordinary word,
whether the token is a word of the gene's approved name (protected: it stays matchable through the full name), and whether the gene is an IntOGen
driver. Classification of candidates is a separate, recorded step (vocab_policy.py, vocab_policy_v5_table.csv).

Output: data/derived/vocab_rescoring/vocab_inventory_screen_candidates.csv
Run: ONCODIGGER_ROOT=<frozen corpus root> python -X utf8 scripts/analysis/vocab_inventory_screen.py
"""
from __future__ import annotations

import csv
import os
import sqlite3
import sys
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ONC = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
WORDS = ROOT / "data" / "reference" / "english_common_words_20k.txt"
WORDS_URL = "https://raw.githubusercontent.com/first20hours/google-10000-english/master/20k.txt"
WORDS_SHA256 = "4ed6e5336d7760d281f7e72df31827da880c861363e820d8c65666b0f10d9ac0"  # the list the paper's screen was run with
INTOGEN = ROOT / "data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"
OUT = ROOT / "data" / "derived" / "vocab_rescoring" / "vocab_inventory_screen_candidates.csv"
ALIAS_MIN, RATIO, SYMBOL_MIN = 200, 5, 300
COLLECTIONS = ["breast", "lung", "colorectal", "prostate", "melanoma", "bladder", "kidney", "pancreas", "liver", "stomach", "esophageal",
               "ovarian", "endometrial", "cervical", "thyroid", "brain", "leukemia", "lymphoma", "myeloma"]


def ensure_words() -> None:
    """Fetch the public word list when it is not on disk, and refuse any copy that differs from the one the paper's screen used."""
    import hashlib
    import urllib.request
    if not WORDS.exists():
        WORDS.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(WORDS_URL, timeout=60) as r:
            WORDS.write_bytes(r.read())
        print(f"fetched {WORDS_URL}", flush=True)
    got = hashlib.sha256(WORDS.read_bytes()).hexdigest()
    if got != WORDS_SHA256:
        sys.exit(f"{WORDS} has SHA-256 {got}, not the recorded {WORDS_SHA256}: the upstream list has changed, so the screen "
                 "would not be the one reported in the paper; stopping")


def counts_for(con, tokens, chunk=400):
    out = {}
    tokens = list(tokens)
    for i in range(0, len(tokens), chunk):
        part = tokens[i:i + chunk]
        q = f"SELECT token, COUNT(DISTINCT pmid) FROM postings WHERE token IN ({','.join('?' * len(part))}) GROUP BY token"
        out.update(con.execute(q, part).fetchall())
    return out


def main() -> None:
    sys.path.insert(0, str(ONC / "src"))
    from oncodigger.search.lexicons import load_lexicons
    single = load_lexicons(ONC)["Genes"]["single"]            # the frozen vocabulary, before any v4/v5 suppression
    phrases = load_lexicons(ONC)["Genes"]["phrases"]
    ensure_words()
    words = set(w.strip().lower() for w in WORDS.read_text(encoding="utf-8").split())
    drivers = {str(s).upper() for s in pd.read_csv(INTOGEN, sep="\t").SYMBOL}
    # approved-name words per gene id (phrases are keyed by anchor token; single-word names appear in `single`)
    name_words = defaultdict(set)
    for anchor, lst in phrases.items():
        for term, e in lst:
            name_words[e["id"]].update(term.split())
    alias_links, symbol_links = [], []
    for tok, ents in single.items():
        for e in ents:
            gid = e["id"]
            if tok == gid.lower():
                if len(tok) >= 3 and tok.isalpha() and tok in words:
                    symbol_links.append((tok, gid))
            elif len(tok) >= 4:
                alias_links.append((tok, gid))
    print(f"alias links {len(alias_links):,} (distinct tokens {len({t for t, _ in alias_links}):,}); ordinary-word symbol tokens {len(symbol_links):,}", flush=True)
    alias_tokens = {t for t, _ in alias_links}
    sym_tokens = {g.lower() for _, g in alias_links} | {t for t, _ in symbol_links}
    best = {}     # (token, gene) -> (alias count, symbol count, collection) for the collection with the largest alias count
    sym_best = {}
    t0 = time.time()
    for coll in COLLECTIONS:
        con = sqlite3.connect(f"file:{(ONC / 'data' / 'pubmed_processed' / f'{coll}_cancer_inclusive_all' / 'lexical_index.db').as_posix()}?mode=ro&immutable=1", uri=True)
        ac = counts_for(con, alias_tokens)
        hot = {t for t, n in ac.items() if n >= ALIAS_MIN}
        sc = counts_for(con, {g.lower() for t, g in alias_links if t in hot} | {t for t, _ in symbol_links})
        for tok, gid in alias_links:
            n = ac.get(tok, 0)
            if n >= ALIAS_MIN:
                ns = sc.get(gid.lower(), 0)
                if n >= RATIO * max(ns, 1) and (tok, gid) not in best or (tok, gid) in best and n > best[(tok, gid)][0] and n >= RATIO * max(ns, 1):
                    best[(tok, gid)] = (n, ns, coll)
        for tok, gid in symbol_links:
            n = sc.get(tok, 0)
            if n >= SYMBOL_MIN and (tok, gid) not in sym_best or (tok, gid) in sym_best and n > sym_best[(tok, gid)][0]:
                if n >= SYMBOL_MIN:
                    sym_best[(tok, gid)] = (n, n, coll)
        con.close()
        print(f"  {coll:11} done ({time.time() - t0:.0f}s): alias candidates so far {len(best)}, symbol candidates {len(sym_best)}", flush=True)
    rows = []
    for (tok, gid), (n, ns, coll) in best.items():
        rows.append({"arm": "alias", "token": tok, "gene": gid, "collection_max": coll, "token_papers": n, "symbol_papers": ns,
                     "ordinary_word": tok in words, "word_of_approved_name": tok in name_words[gid], "intogen_driver": gid in drivers})
    for (tok, gid), (n, ns, coll) in sym_best.items():
        rows.append({"arm": "symbol", "token": tok, "gene": gid, "collection_max": coll, "token_papers": n, "symbol_papers": n,
                     "ordinary_word": True, "word_of_approved_name": tok in name_words[gid], "intogen_driver": gid in drivers})
    rows.sort(key=lambda r: (r["arm"], -r["token_papers"]))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT}: {sum(r['arm'] == 'alias' for r in rows)} alias candidates, {sum(r['arm'] == 'symbol' for r in rows)} symbol candidates", flush=True)


if __name__ == "__main__":
    main()
