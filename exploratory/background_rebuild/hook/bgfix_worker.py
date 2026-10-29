"""bgfix_worker.py -- counts, for one chunk (rowid % chunks == chunk) of one database, the papers in which each entity is found,
using the engine's own per-paper matching block (see bgfix_hook._build_sources). Started by bgfix_hook; not run by hand."""
import importlib
import json
import pickle
import sqlite3
import sys
from collections import Counter
from pathlib import Path

spec = json.load(open(sys.argv[1]))
sys.path[:0] = spec["sys_path"]
sys.path.insert(0, spec["hook_dir"])
import bgfix_hook

module = importlib.import_module(spec["module"])
lookups, patches = pickle.load(open(spec["lookups"], "rb"))
bgfix_hook.apply_patches(module, patches)
passfn, block_sha = bgfix_hook.make_pass_for_worker(module)
con = sqlite3.connect("file:" + Path(spec["db"]).as_posix() + "?mode=ro", uri=True)
K, i = spec["chunks"], spec["chunk"]
seen = set()
n = [0]


def rows():
    for pmid, title, abstract in con.execute("SELECT pmid, title, abstract FROM documents WHERE rowid % ? = ?", (K, i)):
        if pmid in seen:
            continue
        seen.add(pmid)
        n[0] += 1
        if n[0] % 25000 == 0:
            print(i, n[0], flush=True)
        yield {"pmid": pmid}, {"title": title or "", "abstract": abstract or ""}


counts = passfn(rows(), lookups, spec["query_tokens"])
pickle.dump((counts, n[0]), open(spec["out"], "wb"), protocol=4)
print("done", i, n[0], flush=True)
