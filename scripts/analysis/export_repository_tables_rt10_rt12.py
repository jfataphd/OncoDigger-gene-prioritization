"""
export_repository_tables_rt10_rt12.py -- writes Repository Tables RT11 (IntOGen driver recall by cancer at ranks 10/20/50) and RT12
(full-corpus presence breakdown of the 1,740 IntOGen pairs) to data/supplementary_tables/ from the analysis outputs. They were printed
as Supplementary Tables S7 and S10 until the 2026-10-19 slimming of the supplement. RT10 (PubTator disease-entity mapping, formerly
printed Supplementary Table S4) is the printed table itself and is not recomputed.

Inputs : data/derived/query_planes/mutation.csv, the IntOGen reference (intogen_full_rescoring.load_intogen),
         data/derived/intogen_rescoring/intogen_full_corpus_presence_summary.csv
Output : data/supplementary_tables/RT11_*.csv, RT12_*.csv
Run    : python -X utf8 scripts/analysis/export_repository_tables_rt10_rt12.py
"""
import importlib.util
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "supplementary_tables"
spec = importlib.util.spec_from_file_location("ir", ROOT / "scripts" / "analysis" / "intogen_full_rescoring.py")
ir = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ir)


def recall_by_cancer() -> pd.DataFrame:
    _, E = ir.load_intogen()
    E = {(g.upper(), c) for g, c in E}
    rk = pd.read_csv(ROOT / "data" / "derived" / "query_planes" / "mutation.csv")
    rows = []
    for c in ir.STUDY_ORDER:
        drivers = {g for g, cc in E if cc == c}
        top = rk[rk.cancer == c].sort_values("rank").gene.str.upper().tolist()
        row = {"cancer": c, "intogen_driver_pairs_mapped": len(drivers)}
        for k in (10, 20, 50):
            row[f"recall_at_{k}_pct"] = round(100 * len(drivers & set(top[:k])) / len(drivers), 1)
        rows.append(row)
    return pd.DataFrame(rows).sort_values("recall_at_50_pct", kind="stable").reset_index(drop=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    recall_by_cancer().to_csv(OUT / "RT11_IntOGen_recall_by_cancer.csv", index=False)
    p = pd.read_csv(ROOT / "data" / "derived" / "intogen_rescoring" / "intogen_full_corpus_presence_summary.csv")
    p["pct_of_1740"] = p["pct_of_1740"].round(1)
    p.to_csv(OUT / "RT12_full_corpus_presence_breakdown.csv", index=False)
    print("wrote RT11, RT12 to", OUT)


if __name__ == "__main__":
    main()
