"""
proximity_boost_run.py -- post hoc, exploratory proximity-adjusted gene ranking on the frozen engine (Supplementary Results S15).

Retrieval, evidence selection, text handling and the corpus background are those of run_query_planes.py (unchanged). For each gene the
engine credits, the proximity P is the mean over its supporting retrieved papers of 1/(1+d), where d is the number of tokens between the
nearest mention of the gene (symbol, alias or full name in the engine vocabulary) and the nearest query token ('mutation' or 'mutations')
in the same segment (title and abstract are separate segments; adjacent tokens give d = 1; a paper with no segment containing both
contributes 0). The adjusted relevance is  relevance x (1 + w x P)  at full precision (the engine's own 3-decimal rounding of the displayed
Relevance is switched off for this run). w = 0 is the control and reproduces the frozen ranking exactly. The proximity calculation itself is
reference_engine_code/proximity.py, a standalone module copied unchanged from the application checkout (sha256 in
reference_engine_code/README_proximity.md).

Optional argument 'rule' applies the alias-rule sensitivity: single-word aliases credited to two or more genes are removed unless the term is a
gene's own approved symbol, and the PI3K / AKT / PARP family-name aliases are removed.

Env: ONCODIGGER_ROOT = corpus root (verified 2026-09-27 build; for the manuscript run the Lung, Thyroid and Kidney collections are the
MET/REST/FH-corrected copies in data/derived/intogen_rescoring/corrected_corpora).
Usage: PYTHONHASHSEED=0 python -X utf8 scripts/analysis/proximity_boost_run.py <w> <out_dir> [rule]
Outputs: <out_dir>/mutation.csv (cancer, rank, gene; top 50) and <out_dir>/scores_w<w>.csv (cancer, gene, supporting papers,
original relevance, proximity, adjusted relevance for every candidate row).
"""
import builtins
import importlib.util
import os
import sys
from bisect import bisect_left
from pathlib import Path

import pandas as pd

W_LABEL, OUT_DIR = sys.argv[1], sys.argv[2]
W_BOOST = float(W_LABEL)
RULE = "rule" in sys.argv[3:]
ROOT_CORPUS = os.environ["ONCODIGGER_ROOT"]
REPO = Path(__file__).resolve().parents[2]

spec = importlib.util.spec_from_file_location("rqp", REPO / "scripts" / "analysis" / "run_query_planes.py")
m = importlib.util.module_from_spec(spec)
sys.path.insert(0, str(REPO / "scripts" / "analysis"))
spec.loader.exec_module(m)

_sp = importlib.util.spec_from_file_location("proximity_module", REPO / "reference_engine_code" / "proximity.py")
PROX = importlib.util.module_from_spec(_sp)
sys.modules["proximity_module"] = PROX      # the module defines a dataclass and must be registered before it runs
_sp.loader.exec_module(PROX)

QT = {"mutation", "mutations"}
CTX = {}
CANCER_I = [-1]
SCORES = []

if RULE:
    FAMILY = {("pi3k", "PIK3CA"), ("akt", "AKT1"), ("parp", "PARP1")}
    _orig_sup = m.suppress_ambiguous_symbols

    def _sup(lex):
        out = _orig_sup(lex)
        single, dropped = {}, 0
        for t, es in out["single"].items():
            ids = {e["id"] for e in es}
            keep = []
            for e in es:
                if (t, e["id"]) in FAMILY or (len(ids) >= 2 and t != e["id"].lower()):
                    dropped += 1
                    continue
                keep.append(e)
            if keep:
                single[t] = keep
        print("ALIAS RULE: removed", dropped, "term-gene links", flush=True)
        return {"single": single, "phrases": out["phrases"]}

    m.suppress_ambiguous_symbols = _sup

_orig_lm = m.longest_match_lookups


def _lm(lex, discovery):
    res = _orig_lm(lex, discovery)
    CTX["base"] = discovery.full_document_text      # the longest-match text builder; nothing else is changed
    discovery.round = lambda v, n=None: v if n == 3 else builtins.round(v, n)    # full-precision Relevance
    CTX["single"] = res["single"]
    ph = {}
    for a, lst in res["phrases"].items():
        for term, e in lst:
            ph.setdefault(a, []).append((tuple(term.split()), e["id"]))
    CTX["phr"] = ph
    return res


m.longest_match_lookups = _lm


def rank_genes(con, query_text, engine):
    bm25_search, normalize_query, associated_entity_tables, gene_lookups = engine
    tokens = normalize_query(query_text)
    hits = bm25_search(con, tokens, m.SEARCH_POOL, m.OPTIONS)[: m.ANALYSIS_PAPERS]
    tables = associated_entity_tables(con, hits, query_text, (2000, 2026), tokens, gene_lookups,
                                      limit=m.LIMIT, discovery_mode="Primary Evidence", corpus_query="")
    gdf = tables.get("Genes")
    if gdf is None or gdf.empty:
        return []
    cand = []
    for _, row in gdf.iterrows():
        gene = m.GENE_REMAP.get(row["Gene"], row["Gene"])
        cand.append((gene, row["Gene"], float(row["Relevance"]), set(str(row["Supporting PMIDs"]).split(";"))))
    wanted = {c[1] for c in cand} | {c[0] for c in cand}
    pm = sorted({p for c in cand for p in c[3] if p})
    docs = {}
    for i in range(0, len(pm), 500):
        ch = pm[i:i + 500]
        for p, t, a in con.execute(f"SELECT pmid,title,abstract FROM documents WHERE pmid IN ({','.join('?' * len(ch))})", ch):
            docs[str(p)] = (t or "", a or "")
    base, single, phr = CTX["base"], CTX["single"], CTX["phr"]
    terms = {}      # terms of each gene matched in each paper, as the engine's matching loop does (single tokens and phrases)
    for p, (t, a) in docs.items():
        text = base({"title": t, "abstract": a})
        tset = set(text.split())
        padded = " " + text + " "
        for tok in tset:
            for e in single.get(tok, ()):
                if e["id"] in wanted:
                    terms.setdefault((p, e["id"]), set()).add(tok)
            for pt, gid in phr.get(tok, ()):
                if gid in wanted and " " + " ".join(pt) + " " in padded:
                    terms.setdefault((p, gid), set()).add(" ".join(pt))
    CANCER_I[0] += 1
    scored = []
    for gene, raw, rel, sup in cand:
        papers, ebp = [], {}
        for p in sup:
            if p in docs:
                papers.append({"pmid": p, "title": docs[p][0], "abstract": docs[p][1]})
                ebp[p] = terms.get((p, raw), terms.get((p, gene), set()))
        res = PROX.score_supporting_papers(papers, QT, ebp, lambda v: base({"title": v, "abstract": ""}), rel, weight=W_BOOST)
        scored.append((res.adjusted_relevance, gene))
        SCORES.append((m.CANCER_MAP[CANCER_I[0]][1], gene, len(papers), rel, res.proximity, res.adjusted_relevance))
    scored.sort(key=lambda x: -x[0])      # stable: ties keep the engine's order
    ranked, seen = [], set()
    for _, gene in scored:
        if gene in m.GENE_EXCLUDE or gene in seen:
            continue
        try:
            n = con.execute("SELECT COUNT(DISTINCT pmid) FROM postings WHERE token = ?", (gene.lower(),)).fetchone()[0]
        except Exception:
            n = 0
        if n == 0:
            continue
        seen.add(gene)
        ranked.append(gene)
        if len(ranked) >= m.TOP_K:
            break
    return ranked


m.rank_genes = rank_genes
sys.argv = [str(REPO / "scripts" / "analysis" / "run_query_planes.py"), "--planes", "mutation",
            "--oncodigger-root", ROOT_CORPUS, "--out-dir", OUT_DIR, "--allow-nonstandard-build"]
m.main()
pd.DataFrame(SCORES, columns=["cancer", "gene", "supporting_papers", "original_relevance", "proximity", "adjusted_relevance"]).to_csv(
    os.path.join(OUT_DIR, f"scores_w{W_LABEL}.csv"), index=False)
