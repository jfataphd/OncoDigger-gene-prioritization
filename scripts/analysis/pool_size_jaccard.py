"""
pool_size_jaccard.py -- pool-size sensitivity of the top-10 lists (Methods 2.1, Supplementary Methods S1).

The BM25 analysis pool is fixed at the top 1,000 papers per cancer. This script reruns the main pipeline
configuration (same engine setup as intogen_parameter_sensitivity.py) at analysis-pool sizes of 250, 500, 750,
1,000, 1,500 and 2,000 papers (search pool = 1.25 x analysis pool, as in the parameter sweep) for all 19 cancers
and reports the mean top-10 Jaccard similarity (a) of each size with the 1,000-paper pool and (b) between
adjacent sizes.

Output: data/derived/pool_size_jaccard.csv (one row per comparison; also the per-size top-10 lists in
        data/derived/pool_size_top10.csv)
Run:    python -X utf8 scripts/analysis/pool_size_jaccard.py   (needs ONCODIGGER_ROOT; ~25 minutes)
"""
from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import intogen_parameter_sensitivity as ps  # noqa: E402

OUT_DIR = Path(os.environ.get("POOL_SIZE_OUT_DIR", ROOT / "data" / "derived"))
SIZES = [250, 500, 750, 1000, 1500, 2000]


def jaccard(a: list[str], b: list[str]) -> float:
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if a | b else 1.0


def main() -> None:
    if os.environ.get("PYTHONHASHSEED") != "0":
        import subprocess
        sys.exit(subprocess.run([sys.executable, "-X", "utf8"] + sys.argv, env=dict(os.environ, PYTHONHASHSEED="0")).returncode)
    root = ps.DEFAULT_ONCODIGGER
    sys.path.insert(0, str(root / "src"))
    from oncodigger.search import discovery
    from oncodigger.search.discovery import associated_entity_tables
    from oncodigger.search.lexicons import load_lexicons

    lex = load_lexicons(root)
    genes = ps.longest_match_lookups(lex["Genes"], discovery)
    sweep = {"k1": ps.BASELINE_K1, "b": ps.BASELINE_B}
    bm25_search, normalize_query = ps.build_patched_lexical(root, sweep)
    db_base = root / "data" / "pubmed_processed"

    top = {s: {} for s in SIZES}
    for key, label in ps.CANCER_MAP:
        db = next(p for p in (db_base / key / "lexical_index.db", db_base / f"{key}_all" / "lexical_index.db") if p.exists())
        for s in SIZES:
            con = ps.open_corpus(db, discovery)
            top[s][label] = ps.rank_top10(con, bm25_search, normalize_query, associated_entity_tables, genes,
                                          dict(ps.BASELINE_OPTIONS), s, int(s * 1.25))
            con.close()
        print(label, {s: top[s][label][:3] for s in SIZES}, flush=True)

    rows = []
    for s in SIZES:
        if s != 1000:
            rows.append({"comparison": f"{s} vs 1000", "mean_top10_jaccard": float(np.mean([jaccard(top[s][c], top[1000][c]) for c in top[s]]))})
    adj = []
    for a, b in zip(SIZES, SIZES[1:]):
        j = float(np.mean([jaccard(top[a][c], top[b][c]) for c in top[a]]))
        adj.append(j)
        rows.append({"comparison": f"{a}/{b} adjacent", "mean_top10_jaccard": j})
    rows.append({"comparison": "mean of adjacent", "mean_top10_jaccard": float(np.mean(adj))})
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with (OUT_DIR / "pool_size_jaccard.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["comparison", "mean_top10_jaccard"])
        w.writeheader(); w.writerows(rows)
    with (OUT_DIR / "pool_size_top10.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["pool_size", "cancer", "top10"])
        for s in SIZES:
            for c, l in top[s].items():
                w.writerow([s, c, "|".join(l)])
    for r in rows:
        print(f"{r['comparison']}: {r['mean_top10_jaccard']:.3f}")


if __name__ == "__main__":
    main()
