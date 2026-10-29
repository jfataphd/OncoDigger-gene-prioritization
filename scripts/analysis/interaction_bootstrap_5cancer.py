"""
interaction_bootstrap_5cancer.py -- per-cancer breakdown and case-bootstrap CI
for the retrieval x vocabulary interaction, on the matched 5-cancer sample
(Bladder, Esophageal, Endometrial, Myeloma, Thyroid -- the 5 smallest corpora,
chosen for compute-budget reasons before any result was examined, per the
frozen spec amendment log).

Reuses this paper's standing case-bootstrap methodology (resample cancers
with replacement, 2000 replicates, seed 42) rather than inventing new
statistical machinery. Diagnostic only.

Run: python -X utf8 scripts/analysis/interaction_bootstrap_5cancer.py
"""
from __future__ import annotations
import csv
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

INTOGEN_CANCER_MAP = {
    "Bladder Cancer": {"BLCA", "UTUC"}, "Esophageal Cancer": {"ESCA", "ESCC"},
    "Endometrial Cancer": {"UCEC", "UCS"}, "Myeloma": {"PCM"}, "Thyroid Cancer": {"WDTC"},
}
MATCHED_CANCERS = list(INTOGEN_CANCER_MAP.keys())


def load_intogen():
    genes_any, exact_pairs = set(), set()
    path = ROOT / "data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sym = row["SYMBOL"].strip().upper()
            genes_any.add(sym)
            code = row["CANCER_TYPE"].strip()
            for cancer, codes in INTOGEN_CANCER_MAP.items():
                if code in codes:
                    exact_pairs.add((sym, cancer))
    return genes_any, exact_pairs


def top10_from_detail(path, cell_name, cancers):
    rows = list(csv.DictReader(open(path, encoding="utf-8")))
    out = {c: [] for c in cancers}
    for r in rows:
        if r["cell"] != cell_name or r["cancer"] not in out:
            continue
        if int(r["rank"]) <= 10:
            out[r["cancer"]].append((int(r["rank"]), r["gene"]))
    for c in out:
        out[c] = [g for _, g in sorted(out[c])]
    return out


def score(top10_by_cancer, genes_any, exact_pairs):
    any_hits, exact_hits = {}, {}
    for cancer, genes in top10_by_cancer.items():
        any_hits[cancer] = sum(1 for g in genes if g in genes_any)
        exact_hits[cancer] = sum(1 for g in genes if (g, cancer) in exact_pairs)
    return any_hits, exact_hits


def bootstrap_interaction(cell1_any, cell2_any, cell3_any, cell4_any, cancers, n_boot=2000, seed=42):
    rng = random.Random(seed)
    stats = []
    for _ in range(n_boot):
        sample = [rng.choice(cancers) for _ in cancers]
        c1 = sum(cell1_any[c] for c in sample)
        c2 = sum(cell2_any[c] for c in sample)
        c3 = sum(cell3_any[c] for c in sample)
        c4 = sum(cell4_any[c] for c in sample)
        n = len(sample) * 10
        interaction_pp = 100 * (c4 - c2 - c3 + c1) / n
        stats.append(interaction_pp)
    stats.sort()
    lo = stats[int(0.025 * n_boot)]
    hi = stats[int(0.975 * n_boot)]
    point = 100 * (sum(cell4_any[c] for c in cancers) - sum(cell2_any[c] for c in cancers)
                   - sum(cell3_any[c] for c in cancers) + sum(cell1_any[c] for c in cancers)) / (len(cancers) * 10)
    return point, lo, hi


def main():
    genes_any, exact_pairs = load_intogen()

    true_raw_detail = ROOT / "data/derived/intogen_rescoring/true_raw_flat_count_2x2_detail.csv"
    cleaned_detail = ROOT / "data/derived/intogen_rescoring/retrieval_vocab_2x2_faithful_detail.csv"

    cell1_top10 = top10_from_detail(true_raw_detail, "cell1_fullcorpus_true_raw", MATCHED_CANCERS)
    cell3_top10_all = top10_from_detail(true_raw_detail, "cell3_pool_true_raw", MATCHED_CANCERS)
    cell2_top10 = top10_from_detail(cleaned_detail, "cell2_fullcorpus_clean", MATCHED_CANCERS)
    cell4_top10 = top10_from_detail(cleaned_detail, "cell4_pool_clean", MATCHED_CANCERS)

    for name, d in [("cell1 (full corpus, true-raw)", cell1_top10),
                     ("cell2 (full corpus, cleaned)", cell2_top10),
                     ("cell3 (BM25 pool, true-raw)", cell3_top10_all),
                     ("cell4 (BM25 pool, cleaned)", cell4_top10)]:
        missing = [c for c, v in d.items() if len(v) == 0]
        if missing:
            print(f"WARNING: {name} missing data for {missing}")

    cell1_any, cell1_exact = score(cell1_top10, genes_any, exact_pairs)
    cell2_any, cell2_exact = score(cell2_top10, genes_any, exact_pairs)
    cell3_any, cell3_exact = score(cell3_top10_all, genes_any, exact_pairs)
    cell4_any, cell4_exact = score(cell4_top10, genes_any, exact_pairs)

    print("Per-cancer any-cancer hits /10:")
    print(f"{'Cancer':20} {'cell1(raw,full)':>15} {'cell2(clean,full)':>18} {'cell3(raw,pool)':>16} {'cell4(clean,pool)':>18}")
    for c in MATCHED_CANCERS:
        print(f"{c:20} {cell1_any[c]:>15} {cell2_any[c]:>18} {cell3_any[c]:>16} {cell4_any[c]:>18}")
        interaction_c = cell4_any[c] - cell2_any[c] - cell3_any[c] + cell1_any[c]
        print(f"  -> per-cancer any-cancer interaction (cell4-cell2-cell3+cell1): {interaction_c:+d} / 10 slots")

    print()
    point, lo, hi = bootstrap_interaction(cell1_any, cell2_any, cell3_any, cell4_any, MATCHED_CANCERS)
    print(f"ANY-CANCER interaction: {point:+.1f}pp  [95% case-bootstrap CI: {lo:+.1f}, {hi:+.1f}pp]  (n=5 cancers, 2000 reps, seed 42)")

    point_e, lo_e, hi_e = bootstrap_interaction(cell1_exact, cell2_exact, cell3_exact, cell4_exact, MATCHED_CANCERS)
    print(f"EXACT-CONCORDANCE interaction: {point_e:+.1f}pp  [95% case-bootstrap CI: {lo_e:+.1f}, {hi_e:+.1f}pp]")

    print()
    print("Pooled cells (n=50, 5 cancers):")
    for name, any_d, exact_d in [("cell1", cell1_any, cell1_exact), ("cell2", cell2_any, cell2_exact),
                                  ("cell3", cell3_any, cell3_exact), ("cell4", cell4_any, cell4_exact)]:
        a = sum(any_d.values()); e = sum(exact_d.values())
        print(f"  {name}: any={a}/50={100*a/50:.1f}%  exact={e}/50={100*e/50:.1f}%")


if __name__ == "__main__":
    main()
