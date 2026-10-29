"""
make_provenance_record.py -- write data/provenance/outputs_record.csv.

For every derived output (data/derived, data/canonical, figures/png) it records the SHA-256, the size, the
corpus build the outputs of this repository are tied to, and the script(s) in this repository that write a
file of that name. An output with no writer script is labelled "NO GENERATOR IN REPOSITORY" so that
hand-assembled or legacy files are visible instead of silently trusted.

Run after the analysis scripts, then run check_provenance.py:
  python -X utf8 scripts/analysis/make_provenance_record.py
"""
from __future__ import annotations

import csv
import hashlib
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from provenance import file_digest  # noqa: E402
sys_scripts = sorted((ROOT / "scripts").rglob("*.py"))
SOURCES = {p: p.read_text(encoding="utf-8", errors="replace") for p in sys_scripts}
WRITES = re.compile(r"to_csv|\.save\(|savefig|open\(|write_text|json\.dump|DictWriter|csv\.writer|\bsave\(")
BUILD_ID = "2026-09-27"
POLICY_ID = "2026-10-16-v5"   # frozen vocabulary policy v5 applied to every vocabulary-dependent output; corpus-only files are unaffected by it
SKIP_DIRS = {"corrected_corpora", "corrected_corpora_met_strict"}
# .pptx: PowerPoint/python-pptx rewrite the zip container on every save, so the bytes change though no slide content does (verified part by part);
# the figures are recorded through their exported PNGs in figures/png.
SKIP_SUFFIX = {".db", ".pptx"}


# Outputs whose file name is built in code (f-strings, loops) or that are produced by a chain of scripts cannot be linked by name matching.
# Each rule below names the scripts that really write the file; checked against the writing lines of those scripts (2026-10-18).
EXPLICIT = [
    ("data/derived/metastatic_exclusion/secondary_*", ["scripts/analysis/fetch_secondary_pmids.py"]),
    ("data/derived/metastatic_exclusion/build_report_metastasis_excluded*", ["scripts/analysis/build_metastasis_excluded_corpora.py"]),
    ("data/derived/metastatic_exclusion/integrity_check_*", ["scripts/analysis/build_date_truncated_corpora.py"]),                # --verify output
    ("data/derived/metastatic_exclusion/ranking*/mutation.csv", ["scripts/analysis/run_query_planes.py"]),                         # reduced roots, see README
    ("data/derived/metastatic_exclusion/comparison_*", ["scripts/analysis/metastasis_exclusion_summary.py"]),
    ("data/derived/metastatic_exclusion/gene_movements_top10.csv", ["scripts/analysis/metastasis_exclusion_summary.py"]),
    ("data/derived/metastatic_exclusion/liver_selected_genes.csv", ["scripts/analysis/metastasis_exclusion_summary.py"]),
    ("data/derived/metastatic_exclusion/headline_*", ["scripts/analysis/metastasis_exclusion_headline.py"]),
    ("data/derived/rank_aware_metrics.csv", ["scripts/analysis/rank_aware_metrics.py"]),
    ("data/derived/cutoff_matched_rankings/build_report_year_max_*", ["scripts/analysis/build_date_truncated_corpora.py"]),
    ("data/derived/cutoff_matched_rankings/mutation_corpus_to_20??.csv", ["scripts/analysis/run_query_planes.py"]),                 # run on the truncated roots, then copied
    ("data/derived/gene_exclusion_test/mutation_?.csv", ["scripts/analysis/gene_exclusion_readmission_test.py"]),
    ("data/derived/intogen_rescoring/cancermine_intogen_*", ["scripts/analysis/intogen_full_rescoring.py"]),
    ("data/derived/intogen_rescoring/chatgpt_intogen_*", ["scripts/analysis/intogen_full_rescoring.py"]),
    ("data/derived/intogen_rescoring/pooled_control_intogen_*", ["scripts/analysis/intogen_full_rescoring.py"]),
    ("data/derived/intogen_rescoring/diseases_intogen_*", ["scripts/analysis/intogen_diseases_pubtator_rescoring.py"]),
    ("data/derived/intogen_rescoring/pubtator_v2_intogen_*", ["scripts/analysis/intogen_diseases_pubtator_rescoring.py"]),
    ("data/derived/metastatic_exclusion/run_times_batch*", ["hand-written: run() progress files of the two metastasis-exclusion batches (the 2026-09-27 batch driver was not retained)"]),
    ("data/derived/metastatic_exclusion/README.md", ["hand-written documentation of the analysis"]),
]
# Genuinely generator-less files (kept visible on purpose): 
# data/canonical/precision_at_k_detail_pre_vocab_fix.csv.bak (retained backup of an earlier output).


def explicit_producers(rel: str) -> list[str]:
    import fnmatch
    out = []
    for pat, scripts in EXPLICIT:
        if fnmatch.fnmatch(rel, pat):
            out += scripts
    if rel.startswith("figures/png/intogen_") and rel.endswith(".png"):
        stem = rel.rsplit("/", 1)[1][:-4].replace("intogen_", "intogen_build_", 1)
        out += [f"scripts/figures/{stem}.py", "scripts/manuscript/export_pptx_to_png_highres.py"]
    return [x for x in out if x.startswith("exploratory") or x.startswith("hand-written") or (ROOT / x).exists() or not x.startswith("scripts")]


def producers(name: str) -> list[str]:
    stem = name
    out = []
    for p, src in SOURCES.items():
        if p.name == "make_provenance_record.py":
            continue
        if stem in src and WRITES.search(src):
            out.append(str(p.relative_to(ROOT)).replace("\\", "/"))
    return out


def main() -> None:
    rows = []
    files = []
    for base in ("data/derived", "data/canonical", "figures/png", "figures/source"):
        for p in sorted((ROOT / base).rglob("*")):
            if p.is_file() and not (SKIP_DIRS & set(p.parts)) and p.suffix not in SKIP_SUFFIX and p.name != ".DS_Store":
                files.append(p)
    for p in files:
        rel = str(p.relative_to(ROOT)).replace("\\", "/")
        prod = sorted(set(explicit_producers(rel))) or producers(p.name)     # an explicit rule overrides generic name matching (e.g. 'mutation.csv' appears in many readers)
        rows.append({"path": rel, "sha256": file_digest(p), "bytes": p.stat().st_size,
                     "corpus_build": BUILD_ID, "vocabulary_policy": POLICY_ID, "producer_scripts": ";".join(prod) if prod else "NO GENERATOR IN REPOSITORY"})
    dst = ROOT / "data" / "provenance" / "outputs_record.csv"
    dst.parent.mkdir(parents=True, exist_ok=True)
    with dst.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["path", "sha256", "bytes", "corpus_build", "vocabulary_policy", "producer_scripts"])
        w.writeheader()
        w.writerows(rows)
    nogen = [r["path"] for r in rows if r["producer_scripts"].startswith("NO GENERATOR")]
    print(f"Recorded {len(rows)} outputs; {len(nogen)} have no generator script in this repository:")
    for x in nogen:
        print("  ", x)
    print("Wrote", dst)


if __name__ == "__main__":
    main()
