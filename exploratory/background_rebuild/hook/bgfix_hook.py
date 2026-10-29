"""
bgfix_hook.py -- run any UNCHANGED analysis script with the exact background count (docs/BACKGROUND_COUNT_SPEC.md).

How it works
  * Activated by sitecustomize.py when BGFIX_ENABLE=1 (nothing else in the script or the engine file is edited).
  * When a module named *.discovery is imported, its associated_entity_tables() source is read and exactly two things are changed:
      1. one line is inserted after `documents = _docs_for_results(index, results)` that starts the exact-background context;
      2. the line `background_count = entity_background_papers(index, list(matched_terms[entity_id])) or support_count`
         is replaced by a lookup of the exact count.
    Everything else (matching, weights, formula, multipliers, ordering) is the engine's own code. If either line is not found
    exactly once, the hook raises; it never silently falls back to the old count.
  * The exact count B for an entity is the number of distinct papers of the collection in which THE ENGINE'S OWN per-paper matching
    block (extracted from the same function source) finds the entity, with the lookups that the calling script passed. So S and B use
    one recogniser, whichever vocabulary the script uses (cleaned, true-raw, unfiltered, ...).
  * Counts are built once per (database content, lookups content, query, matching code) by BGFIX_WORKERS separate worker processes
    (bgfix_worker.py) and stored with that stamp; a stored count is refused if any part of the stamp differs.
  * Guards on every call: the pool support recomputed with the extracted block must equal the engine's support_count for every entity;
    S <= B <= collection size; 0 < E <= C/N.
Environment: BGFIX_ENABLE=1, BGFIX_CACHE=<dir>, BGFIX_WORKERS=<n, default 14>, BGFIX_STATS=<json file written at exit>.
"""
import atexit
import hashlib
import importlib.abc
import importlib.machinery
import inspect
import json
import marshal
import os
import pickle
import re
import types
import subprocess
import sys
import tempfile
import textwrap
import time
from collections import Counter, defaultdict
from pathlib import Path

HOOK_VERSION = "bgfix-hook-1"
HERE = Path(__file__).resolve().parent
STATS = {"hook_version": HOOK_VERSION, "calls": 0, "entities": 0, "engine_b_lower": 0, "engine_b_equal": 0, "engine_b_higher": 0,
         "floor_would_bind": 0, "tables": {}}

_START = "documents = _docs_for_results(index, results)"
_LOOP = "for rank, result in enumerate(results):"
_END = "enriched_entities = []"
_BLINE = "background_count = entity_background_papers(index, list(matched_terms[entity_id])) or support_count"


def _once(text, needle):
    if text.count(needle) != 1:
        raise RuntimeError(f"bgfix: expected exactly one occurrence of {needle!r}, found {text.count(needle)}")


def _build_sources(fn_src: str):
    """Return (patched associated_entity_tables source, pass-function source) built from the engine's own function text."""
    for n in (_START, _LOOP, _END, _BLINE):
        _once(fn_src, n)
    lines = fn_src.splitlines()
    i_start = next(i for i, l in enumerate(lines) if l.strip() == _START)
    i_loop = next(i for i, l in enumerate(lines) if l.strip() == _LOOP)
    i_end = next(i for i, l in enumerate(lines) if l.strip() == _END)
    i_q = next(i for i, l in enumerate(lines) if l.strip().startswith("query_terms = {"))
    ind = lambda l: len(l) - len(l.lstrip())
    # ---- patched engine function
    pl = list(lines)
    pl[i_start] = lines[i_start] + "\n" + " " * ind(lines[i_start]) + "_bgfix_ctx = _bgfix_begin(index, lookups, query_tokens, results, documents, corpus_doc_count_fn=lambda: _doc_count(index))"
    k = next(i for i, l in enumerate(pl) if l.strip() == _BLINE)
    pl[k] = " " * ind(pl[k]) + "background_count = _bgfix_b(_bgfix_ctx, entity_id, labels[entity_id], support_count, len(results), corpus_doc_count, entity_background_papers(index, list(matched_terms[entity_id])) or support_count)"
    patched = "\n".join(pl)
    # ---- pass function: the engine's own per-paper block, counting papers per entity
    head = lines[i_q:i_start]
    loop = lines[i_loop:i_end]
    out = []
    for j, l in enumerate(loop):
        s = l.strip()
        if j == 0:
            out.append(" " * ind(l) + "for result, doc in rows:")
        elif s.startswith("doc = documents.get(result["):
            continue
        elif s.startswith("weight = 1.0 / math.sqrt(rank + 1)"):
            continue
        elif s == "scores[entity_id] += weight":
            out.append(" " * ind(l) + "counts[entity_id] += 1")
        elif s == 'evidence[entity_id].add(result["pmid"])':
            out.append(" " * ind(l) + "pass")
        else:
            out.append(l)
    body = "\n".join(out)
    if re.search(r"\brank\b|\bweight\b", body):
        raise RuntimeError("bgfix: engine block uses rank/weight outside the expected lines; extraction not safe")
    if "counts[entity_id] += 1" not in body:
        raise RuntimeError("bgfix: scoring line not found in engine block")
    pass_src = ("def _bgfix_pass(rows, lookups, query_tokens):\n"
                + "\n".join(head) + "\n    counts = Counter()\n" + body + "\n    return counts\n")
    return patched, pass_src, hashlib.sha256(body.encode()).hexdigest()


def _sha256_file(path):
    p = Path(path)
    st = p.stat()
    reg = Path(os.environ["BGFIX_CACHE"]) / "_dbhash.json"
    try:
        d = json.load(open(reg))
    except Exception:
        d = {}
    key = str(p.resolve())
    if key in d and d[key][0] == st.st_size and d[key][1] == st.st_mtime_ns:
        return d[key][2]
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while True:
            b = fh.read(1 << 24)
            if not b:
                break
            h.update(b)
    d[key] = [st.st_size, st.st_mtime_ns, h.hexdigest()]
    fd, tmp = tempfile.mkstemp(dir=reg.parent, suffix=".tmp")
    with os.fdopen(fd, "w") as fh:
        json.dump(d, fh)
    os.replace(tmp, reg)
    return h.hexdigest()


def lookups_hash(lookups) -> str:
    h = hashlib.sha256()
    for bucket in sorted(lookups):
        lk = lookups[bucket]
        h.update(repr(("bucket", bucket)).encode())
        for token, ents in sorted(lk.get("single", {}).items()):
            h.update(repr((token, sorted((str(e["id"]), e["label"], e.get("category", "")) for e in ents))).encode())
        for anchor, pes in sorted(lk.get("phrases", {}).items()):
            h.update(repr((anchor, sorted((p, str(e["id"]), e["label"], e.get("category", "")) for p, e in pes))).encode())
    return h.hexdigest()


def _is_pristine(v, pristine):
    return any(v is o for o in pristine.values())


def _ser_func(f, pristine):
    if f.__defaults__ or f.__kwdefaults__:
        raise RuntimeError(f"bgfix: patched function {f.__name__} has defaults; cannot be replayed in workers")
    free = {}
    for name, cell in zip(f.__code__.co_freevars, f.__closure__ or ()):
        v = cell.cell_contents
        if inspect.isfunction(v) and _is_pristine(v, pristine):
            free[name] = ("pristine", next(k for k, o in pristine.items() if o is v))
        elif inspect.isfunction(v):
            free[name] = ("func", _ser_func(v, pristine))
        elif isinstance(v, re.Pattern):
            free[name] = ("re", v.pattern, v.flags)
        else:
            try:
                free[name] = ("val", pickle.dumps(v, protocol=4))
            except Exception as exc:
                raise RuntimeError(f"bgfix: cannot serialise closure value {name} of {f.__name__}: {exc}")
    return {"name": f.__name__, "code": marshal.dumps(f.__code__), "free": free}


def serialise_patches(module):
    """Engine-module functions that the calling script replaced after import (for example the longest-match masking of
    full_document_text). Workers must run the same functions or they would match differently from the parent."""
    pristine = _STATE[module.__name__]["pristine"]
    out = {}
    for name, orig in pristine.items():
        cur = module.__dict__.get(name)
        if cur is not orig:
            if not inspect.isfunction(cur):
                raise RuntimeError(f"bgfix: {name} replaced by a non-function")
            out[name] = _ser_func(cur, pristine)
    return out


def _build_func(spec, module, pristine):
    cells = []
    for name in marshal.loads(spec["code"]).co_freevars:
        kind, *val = spec["free"][name]
        if kind == "pristine":
            v = pristine[val[0]]
        elif kind == "func":
            v = _build_func(val[0], module, pristine)
        elif kind == "re":
            v = re.compile(val[0], val[1])
        else:
            v = pickle.loads(val[0])
        cells.append(types.CellType(v))
    return types.FunctionType(marshal.loads(spec["code"]), module.__dict__, spec["name"], None, tuple(cells) or None)


def apply_patches(module, patches):
    pristine = {n: o for n, o in vars(module).items() if inspect.isfunction(o) and o.__module__ == module.__name__
                and n != "associated_entity_tables" and not n.startswith("_bgfix")}
    built = {name: _build_func(spec, module, pristine) for name, spec in patches.items()}
    module.__dict__.update(built)


def pristine_functions(module):
    return {n: o for n, o in vars(module).items() if inspect.isfunction(o) and o.__module__ == module.__name__
            and n != "associated_entity_tables" and not n.startswith("_bgfix")}


class _Ctx:
    pass


_STATE = {}   # per-module: pass function; per (db, lookups id): table


def _get_pass(module):
    return _STATE[module.__name__]["pass"]


def _db_path(index):
    return next(r[2] for r in index.execute("PRAGMA database_list").fetchall() if r[1] == "main")


def _table(module, index, lookups, query_tokens, n_docs):
    """Exact B per (bucket, entity id) for this database + lookups + query; built once, stored with a stamp."""
    db = _db_path(index)
    memo_key = (db, id(lookups), tuple(query_tokens))
    st = _STATE[module.__name__]
    if memo_key in st["memo"]:
        return st["memo"][memo_key]
    patches = serialise_patches(module)
    cache = Path(os.environ["BGFIX_CACHE"])
    cache.mkdir(parents=True, exist_ok=True)
    stamp = {"hook_version": HOOK_VERSION, "db_sha256": _sha256_file(db), "db_name": Path(db).parent.name,
             "lookups_sha256": lookups_hash(lookups), "module_patches_sha256": hashlib.sha256(pickle.dumps(patches, protocol=4)).hexdigest(), "query_tokens": list(query_tokens), "engine_block_sha256": st["block_sha"],
             "doc_count": n_docs}
    name = hashlib.sha1(json.dumps(stamp, sort_keys=True).encode()).hexdigest()[:20]
    fcsv, fst = cache / f"{name}.json", cache / f"{name}.stamp.json"
    if fst.exists() and fcsv.exists():
        got = json.load(open(fst))
        if got != stamp:
            raise RuntimeError(f"bgfix: stale count table {name}")
        raw = json.load(open(fcsv))
        table = {(b, i): v for b, i, v in raw["counts"]}
        n = raw["papers"]
    else:
        n, table = _build(module, db, lookups, query_tokens, name, n_docs, patches)
        fd, tmp = tempfile.mkstemp(dir=cache, suffix=".tmp")
        with os.fdopen(fd, "w") as fh:
            json.dump({"papers": n, "counts": [[b, i, v] for (b, i), v in sorted(table.items(), key=lambda kv: (str(kv[0][0]), str(kv[0][1])))]}, fh)
        os.replace(tmp, fcsv)
        fd, tmp = tempfile.mkstemp(dir=cache, suffix=".tmp")
        with os.fdopen(fd, "w") as fh:
            json.dump(stamp, fh, indent=1, sort_keys=True)
        os.replace(tmp, fst)
    if n != n_docs:
        raise RuntimeError(f"bgfix: counted {n} papers but the engine's corpus size is {n_docs}")
    STATS["tables"][name] = {"db": stamp["db_name"], "papers": n, "entities": len(table)}
    st["memo"][memo_key] = (table, n)
    return st["memo"][memo_key]


def _build(module, db, lookups, query_tokens, name, n_docs, patches):
    K = int(os.environ.get("BGFIX_WORKERS", "14"))
    cache = Path(os.environ["BGFIX_CACHE"])
    work = cache / f"_work_{name}"
    work.mkdir(parents=True, exist_ok=True)
    lk = work / "lookups.pkl"
    with open(lk, "wb") as fh:
        pickle.dump((lookups, patches), fh, protocol=4)
    procs = []
    t0 = time.time()
    for i in range(K):
        spec = {"sys_path": [p for p in sys.path if p], "module": module.__name__, "db": db, "lookups": str(lk),
                "query_tokens": list(query_tokens), "chunk": i, "chunks": K, "out": str(work / f"part{i}.pkl"),
                "hook_dir": str(HERE)}
        sp = work / f"spec{i}.json"
        json.dump(spec, open(sp, "w"))
        env = {k: v for k, v in os.environ.items() if k != "BGFIX_ENABLE"}
        env["PYTHONHASHSEED"] = "0"
        env["PYTHONPATH"] = ""
        log = open(work / f"worker{i}.log", "w")
        procs.append((subprocess.Popen([sys.executable, "-X", "utf8", str(HERE / "bgfix_worker.py"), str(sp)], env=env, stdout=log, stderr=subprocess.STDOUT), log))
    codes = []
    for p, log in procs:
        codes.append(p.wait())
        log.close()
    if any(codes):
        raise RuntimeError(f"bgfix: worker failure {codes}; see {work}")
    total, n = Counter(), 0
    for i in range(K):
        with open(work / f"part{i}.pkl", "rb") as fh:
            c, k = pickle.load(fh)
        total.update(c)
        n += k
    print(f"[bgfix] built {name} for {Path(db).parent.name}: {n} papers, {len(total)} entities, {time.time() - t0:.0f}s", file=sys.stderr, flush=True)
    for f in work.iterdir():
        f.unlink()
    work.rmdir()
    return n, dict(total)


def _begin(module):
    def begin(index, lookups, query_tokens, results, documents, corpus_doc_count_fn):
        n_docs = corpus_doc_count_fn()
        table, n = _table(module, index, lookups, query_tokens, n_docs)
        rows = [(r, documents.get(r["pmid"], {})) for r in results]
        pool = _get_pass(module)(((r, {"title": d.get("title", ""), "abstract": d.get("abstract", "")}) for r, d in rows), lookups, query_tokens)
        c = _Ctx()
        c.table, c.n, c.pool = table, n, pool
        STATS["calls"] += 1
        return c
    return begin


def _b(ctx, entity_id, label, support, n_pool, corpus_n, engine_b):
    if ctx.pool.get(entity_id, 0) != support:
        raise RuntimeError(f"bgfix STOP: pool support differs for {entity_id} {label}: engine {support}, extracted block {ctx.pool.get(entity_id, 0)}")
    b = ctx.table.get(entity_id)
    if b is None:
        raise RuntimeError(f"bgfix STOP: no exact count for {entity_id} {label}")
    if not (support <= b <= ctx.n) or ctx.n != corpus_n:
        raise RuntimeError(f"bgfix STOP: invariant violated for {label}: S={support} B={b} C={ctx.n} engineC={corpus_n}")
    E = (support / n_pool) / (b / corpus_n)
    if not (0 < E <= corpus_n / n_pool + 1e-9):
        raise RuntimeError(f"bgfix STOP: enrichment {E} outside (0, C/N] for {label}")
    STATS["entities"] += 1
    if engine_b < b:
        STATS["engine_b_lower"] += 1
    elif engine_b == b:
        STATS["engine_b_equal"] += 1
    else:
        STATS["engine_b_higher"] += 1
    if engine_b < support:
        STATS["floor_would_bind"] += 1
    return b


def patch_module(module):
    if getattr(module, "_bgfix_patched", False) or not hasattr(module, "associated_entity_tables"):
        return
    fn_src = textwrap.dedent(inspect.getsource(module.associated_entity_tables))
    patched, pass_src, block_sha = _build_sources(fn_src)
    pristine = pristine_functions(module)
    ns = module.__dict__
    ns["_bgfix_begin"] = _begin(module)
    ns["_bgfix_b"] = _b
    exec(compile(pass_src, f"<bgfix pass {module.__name__}>", "exec"), ns)
    _STATE[module.__name__] = {"pass": ns["_bgfix_pass"], "block_sha": block_sha, "memo": {}, "pristine": pristine}
    exec(compile(patched, f"<bgfix patched {module.__name__}>", "exec"), ns)
    module._bgfix_patched = True
    print(f"[bgfix] patched {module.__name__} (engine block {block_sha[:12]})", file=sys.stderr, flush=True)


def make_pass_for_worker(module):
    fn_src = textwrap.dedent(inspect.getsource(module.associated_entity_tables))
    _, pass_src, block_sha = _build_sources(fn_src)
    ns = module.__dict__
    exec(compile(pass_src, f"<bgfix pass {module.__name__}>", "exec"), ns)
    return ns["_bgfix_pass"], block_sha


class _Finder(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path, target=None):
        if not fullname.endswith(".discovery") or fullname in _Finder.busy:
            return None
        _Finder.busy.add(fullname)
        try:
            spec = importlib.machinery.PathFinder.find_spec(fullname, path)
        finally:
            _Finder.busy.discard(fullname)
        if spec is None or spec.loader is None:
            return None
        orig = spec.loader.exec_module

        def exec_module(module, _orig=orig):
            _orig(module)
            patch_module(module)
        spec.loader.exec_module = exec_module
        return spec


_Finder.busy = set()


def install():
    sys.meta_path.insert(0, _Finder())

    def dump():
        p = os.environ.get("BGFIX_STATS")
        if p:
            json.dump(STATS, open(p, "w"), indent=1)
    atexit.register(dump)
