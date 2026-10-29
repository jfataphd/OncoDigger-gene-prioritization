"""
title_case_blind_spot_check.py -- how often the capitalization filter retains a posting on a mixed-case form
only (Scope and limitations, "Capitalization filter and article-title conventions").

The two-step false-positive rule (vocab_collision_pipeline.capitalization_and_negative_filter) treats any
standalone, non-sentence-initial occurrence of a symbol that is not all-lowercase as evidence of genuine gene
usage. Journal titles are often Title Cased, so an ordinary word ("Met", "Rest") capitalized only by
title-casing convention can pass. This script bounds that blind spot for the two symbols whose ordinary-word
senses matter most, MET and REST, across all 19 corpora.

For each corpus and symbol it takes the postings that the rule RETAINS (capitalization step passed, no
negative phrase), then counts how many of those retained postings (a) contain NO all-caps occurrence of the symbol (their only qualifying matches are
mixed-case forms such as "Met" or "Rest"), and, more tightly, (b) have every qualifying match inside a
Title-Cased article title (at least 70% of the title's words of four or more letters start with a capital
letter). Column (b) counts the postings that a title-casing convention alone could have admitted. Both are upper
bounds on the effect, not counts of confirmed errors: a mixed-case match can also be a genuine gene mention.

Output: data/derived/vocab_rescoring/title_case_blind_spot_check.csv
Run: python -X utf8 scripts/analysis/title_case_blind_spot_check.py   (needs ONCODIGGER_ROOT)
"""
from __future__ import annotations

import csv
import json
import re
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from vocab_collision_pipeline import NEGATIVE_PHRASES, CANCER_KEY, BASE_DB, BASE_RAW  # noqa: E402

OUT = ROOT / "data" / "derived" / "vocab_rescoring" / "title_case_blind_spot_check.csv"
SENTENCE_START = re.compile(r"(^|[.!?]\s+)\s*$")
SYMBOLS = ["MET", "REST"]


def main() -> None:
    rows = []
    for gene in SYMBOLS:
        token = gene.lower()
        pattern = re.compile(r"\b" + re.escape(token) + r"\b", re.IGNORECASE)
        neg = NEGATIVE_PHRASES.get(gene, [])
        for cancer, key in CANCER_KEY.items():
            con = sqlite3.connect(f"file:{(BASE_DB / (key + '_all') / 'lexical_index.db').as_posix()}?mode=ro", uri=True)
            pmids = {str(p[0]) for p in con.execute("SELECT DISTINCT pmid FROM postings WHERE token = ?", (token,))}
            con.close()
            retained = mixed_only = title_only = 0
            with open(BASE_RAW / key / "abstracts.jsonl", encoding="utf-8") as f:
                for line in f:
                    rec = json.loads(line)
                    pmid = str(rec.get("pmid", ""))
                    if pmid not in pmids:
                        continue
                    text = (rec.get("title") or "") + " " + (rec.get("abstract") or "")
                    title = rec.get("title") or ""
                    quals = [m.group(0) for m in pattern.finditer(text)
                             if m.group(0) != token and not SENTENCE_START.search(text[:m.start()])]
                    in_title = [m for m in pattern.finditer(title)
                                if m.group(0) != token and not SENTENCE_START.search(title[:m.start()])]
                    if not quals:
                        continue
                    low = text.lower()
                    if neg and any(p in low for p in neg):
                        continue
                    retained += 1
                    if not any(q == gene for q in quals):
                        mixed_only += 1
                        words = [w for w in re.findall(r"[A-Za-z]+", title) if len(w) >= 4]
                        titlecased = bool(words) and sum(w[0].isupper() for w in words) / len(words) >= 0.7
                        if titlecased and len(in_title) == len(quals):
                            title_only += 1
            rows.append({"symbol": gene, "cancer": cancer, "postings": len(pmids), "retained": retained,
                         "retained_mixed_case_only": mixed_only,
                         "retained_title_cased_title_only": title_only})
    with OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for gene in SYMBOLS:
        r = [x for x in rows if x["symbol"] == gene]
        ret, mix = sum(x["retained"] for x in r), sum(x["retained_mixed_case_only"] for x in r)
        tit = sum(x["retained_title_cased_title_only"] for x in r)
        print(f"{gene}: {ret} retained postings across 19 corpora; {mix} ({100 * mix / max(ret, 1):.1f}%) contain only mixed-case forms; "
              f"{tit} ({100 * tit / max(ret, 1):.1f}%) are admitted only by a mixed-case match in a Title-Cased title")
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
