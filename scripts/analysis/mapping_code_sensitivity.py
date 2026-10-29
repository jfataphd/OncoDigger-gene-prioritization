"""
mapping_code_sensitivity.py -- does the exact cancer-type concordance (142/190) depend on the broader IntOGen codes in the mapping?

Recomputes the OncoDigger top-10 exact and any-cancer concordance with IntOGen cancer codes removed from the mapping of Supplementary
Methods S5, one code at a time and in groups, and lists the gene-cancer pairs whose exact match rests only on a removed code.
Uses the same loader and mapping as intogen_full_rescoring.py. No engine run.
Output: data/derived/mapping_code_sensitivity.csv
Run   : python -X utf8 scripts/analysis/mapping_code_sensitivity.py
"""
import csv
import importlib.util
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("ir", ROOT / "scripts" / "analysis" / "intogen_full_rescoring.py")
ir = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ir)

top10 = ir.load_top10("oncodigger", ir.SYSTEMS["oncodigger"])
rows = list(csv.DictReader(open(ir.INTOGEN_DRIVERS, encoding="utf-8"), delimiter="\t"))


def pairs_without(drop):
    p = set()
    for r in rows:
        if r["CANCER_TYPE"] in drop:
            continue
        for cancer, codes in ir.CANCER_MAP.items():
            if r["CANCER_TYPE"] in codes:
                p.add((r["SYMBOL"], cancer))
    return p


def exact(drop):
    p = pairs_without(drop)
    return int(sum((g, c) in p for g, c in zip(top10.gene, top10.cancer)))


base = exact(set())
print("baseline exact", base)
out = []
tests = [("UTUC",), ("MDS",), ("CHOL",), ("UTUC", "MDS", "CHOL")]
# the other codes that are not the main histology of their study cancer
tests += [(c,) for c in ("UM", "WT", "PANET", "UCS", "CEAD", "CML", "PAST", "MBL", "EPM", "ATRT", "LIHB", "SCLC")]
tests.append(("UTUC", "MDS", "CHOL", "UM", "WT", "PANET", "UCS", "CEAD", "CML", "PAST", "MBL", "EPM", "ATRT", "LIHB", "SCLC"))
full = pairs_without(set())
for drop in tests:
    p = pairs_without(set(drop))
    lost = [(g, c) for g, c in zip(top10.gene, top10.cancer) if (g, c) in full and (g, c) not in p]
    out.append({"codes_removed": "+".join(drop), "exact": exact(set(drop)), "lost_pairs": "; ".join(f"{g} ({c})" for g, c in lost)})
res = pd.DataFrame(out)
res.to_csv(ROOT / "data" / "derived" / "mapping_code_sensitivity.csv", index=False)
print(res.to_string(index=False))
