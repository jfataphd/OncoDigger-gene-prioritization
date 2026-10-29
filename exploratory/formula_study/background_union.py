"""
background_union.py -- background document counts that include every active alias of a gene (exploratory).

The engine's background count B is the smallest postings count among the gene's matched single-word terms and is 0 for genes matched only as
a phrase (for example NKX2-1, normalised 'nkx2 1'), in which case it falls back to the gene's own supporting count. Here B_union is the number of
collection papers that contain ANY active single-word alias of the gene in the postings (phrase-only genes: normalised-text phrase match).
Reads features_raw.csv; needs the engine and the corrected databases. Output: background_union.csv (cancer, gene, B_union, basis)
Run: PYTHONHASHSEED=0 python -X utf8 background_union.py <ONCODIGGER_ROOT>
"""
import importlib.util, os, re, sys
from pathlib import Path
import pandas as pd
assert os.environ.get("PYTHONHASHSEED") == "0"
HERE = Path(__file__).resolve().parent; REPO = HERE.parents[1]; ROOT = Path(sys.argv[1])
sys.path.insert(0, str(REPO / "scripts" / "analysis"))
spec = importlib.util.spec_from_file_location("rqp", REPO / "scripts" / "analysis" / "run_query_planes.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
sys.path.insert(0, str(ROOT / "src"))
from oncodigger.search import discovery
from oncodigger.search.lexicons import load_lexicons, norm_text
lex = load_lexicons(ROOT)
G = m.longest_match_lookups(m.suppress_ambiguous_symbols(lex["Genes"]), discovery)
single, phrase = {}, {}
for tok, ents in G["single"].items():
    for e in ents:
        single.setdefault(e["id"], set()).add(tok)
for anchor, lst in G["phrases"].items():
    for ph, e in lst:
        phrase.setdefault(e["id"], set()).add(ph)
f = pd.read_csv(HERE / "features_raw.csv")
inv = {v: k for k, v in m.GENE_REMAP.items()}
rows = []
for key, label in m.CANCER_MAP:
    con = m.open_corpus(ROOT / "data" / "pubmed_processed" / f"{key}_all" / "lexical_index.db", discovery)
    genes = f[f.cancer == label].gene.tolist()
    need_text = []
    for g in genes:
        eids = {g, inv.get(g, g)}
        toks = set().union(*[single.get(e, set()) for e in eids])
        if toks:
            q = ",".join("?" * len(toks))
            n = con.execute(f"SELECT COUNT(DISTINCT pmid) FROM postings WHERE token IN ({q})", list(toks)).fetchone()[0]
            rows.append({"cancer": label, "gene": g, "B_union": n, "basis": "single-word aliases"})
        else:
            need_text.append((g, set().union(*[phrase.get(e, set()) for e in eids])))
    if need_text:
        texts = [" " + norm_text(f"{t or ''} {a or ''}") + " " for t, a in con.execute("SELECT title, abstract FROM documents")]
        for g, phs in need_text:
            n = sum(1 for t in texts if any(f" {p} " in t for p in phs)) if phs else 0
            rows.append({"cancer": label, "gene": g, "B_union": n, "basis": "phrase text scan"})
    con.close(); print(label, len(genes), len(need_text), flush=True)
pd.DataFrame(rows).to_csv(HERE / "background_union.csv", index=False)
print("done")
