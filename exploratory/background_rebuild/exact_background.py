"""
exact_background.py -- exact background counts for one or more collections (see docs/BACKGROUND_COUNT_SPEC.md).
Writes b_exact_<name>.csv + b_exact_<name>.stamp.json (atomic, version-stamped). Stops on a violated invariant.
Run: PYTHONHASHSEED=0 python -X utf8 exact_background.py <ONCODIGGER_ROOT> "Lung Cancer" ...      (pooled index: name "POOLED")
"""
import os
import sys
import time
from pathlib import Path

assert os.environ.get("PYTHONHASHSEED") == "0", "set PYTHONHASHSEED=0"
sys.path.insert(0, str(Path(__file__).resolve().parent))
import bg_matcher as bg

HERE = Path(__file__).resolve().parent
ROOT = Path(sys.argv[1])
mt = bg.Matcher(ROOT)
m, d = mt.m, mt.discovery
from oncodigger.search.lexical import bm25_search, normalize_query
tokens = normalize_query(m.PLANES["mutation"])
targets = [(k, l) for k, l in m.CANCER_MAP if l in sys.argv[2:]]
if "POOLED" in sys.argv[2:]:
    targets.append(("pooled_common_cancers", "POOLED"))
for key, label in targets:
    t0 = time.time()
    dbp = ROOT / "data" / "pubmed_processed" / f"{key}_all" / "lexical_index.db"
    st = bg.stamp(mt, dbp)
    con = m.open_corpus(dbp, d)
    hits = bm25_search(con, tokens, m.SEARCH_POOL, m.OPTIONS)[: m.ANALYSIS_PAPERS] if label != "POOLED" else []
    pool_ids = {str(h["pmid"]) for h in hits}
    b, s, n = mt.count(con.execute("SELECT pmid, title, abstract FROM documents"), pool_ids,
                       progress=lambda k: print(label, k, f"{time.time() - t0:.0f}s", flush=True))
    con.close()
    bad = [(k, s[k], b[k]) for k in s if s[k] > b[k]]
    if bad or any(v > n for v in b.values()):
        print("STOP: invariant violated", label, bad[:5], flush=True)
        sys.exit(2)
    st["papers"] = n
    rows = [[eid, lab, v, s.get((eid, lab), 0), n] for (eid, lab), v in sorted(b.items())]
    bg.save(HERE / f"b_exact_{label.replace(' ', '_')}", ["entity_id", "label", "B_exact", "S_pool", "collection_papers"], rows, st)
    print(label, "papers", n, "entities", len(b), f"{time.time() - t0:.0f}s", flush=True)
