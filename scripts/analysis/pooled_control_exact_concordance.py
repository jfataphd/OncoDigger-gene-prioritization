"""
pooled_control_exact_concordance.py — the (gene, cancer) pairs of the cancer-blind pooled control.

The pooled control (all 19 collections pooled, one ranking) returns the same top 10 genes for every cancer. This lists those
10 genes once for each of the 19 cancers, 190 (pooled gene x cancer) pairs, so that the control is scored against each cancer
exactly as the cancer-specific rankings are. The file is read by intogen_full_rescoring.py (IntOGen any-cancer and exact
concordance of the pooled control) and by verify_headline_results.py (corpus-separation effect). The file name is historical;
it holds no concordance values.

Run: python -X utf8 scripts/analysis/pooled_control_exact_concordance.py
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

CANCERS = [
    "Bladder Cancer", "Brain Cancer", "Breast Cancer", "Cervical Cancer", "Colorectal Cancer",
    "Endometrial Cancer", "Esophageal Cancer", "Kidney Cancer", "Leukemia", "Liver Cancer",
    "Lung Cancer", "Lymphoma", "Melanoma", "Myeloma", "Ovarian Cancer", "Pancreatic Cancer",
    "Prostate Cancer", "Stomach Cancer", "Thyroid Cancer",
]


def main() -> None:
    pooled = pd.read_csv(ROOT / "data" / "derived" / "pooled_control_mutation.csv")
    pooled_top10 = pooled[pooled["rank"] <= 10]["gene"].tolist()

    rows = [{"cancer": c, "gene": g} for c in CANCERS for g in pooled_top10]
    df = pd.DataFrame(rows)
    print(f"Pooled top-10 genes: {pooled_top10}")
    print(f"Total (gene, cancer) pairs: {len(df)}")
    out = ROOT / "data" / "derived" / "pooled_control_exact_concordance.csv"
    df.to_csv(out, index=False)
    print(f"Saved {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
