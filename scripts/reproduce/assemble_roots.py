"""
assemble_roots.py -- build the corpus "roots" that the analysis scripts expect, from the 19 corpus databases and files in this repository.

A root is a folder holding data/pubmed_processed/<collection>_all/lexical_index.db plus the four lexicon CSVs and src/oncodigger/search (the
frozen engine). Nothing here is downloaded; the engine and lexicons are the vendored copies in reference_engine_code/.

  frozen   ONCODIGGER_ROOT of the main analysis: 19 collections linked from --corpora.
  corr     as frozen, except Lung, Thyroid and Kidney point at the MET / REST / FH-corrected copies
           (data/derived/intogen_rescoring/corrected_corpora, built by vocab_corrected_rerun.py). Used by the proximity runs.

Usage: python -X utf8 scripts/reproduce/assemble_roots.py --engine-into <root> [<root> ...]   (lexicons + engine only)
       python -X utf8 scripts/reproduce/assemble_roots.py --corpora <dir with the 19 *_cancer_inclusive_all folders> --work <work dir> [--raw <dir with pubmed_raw>] [--corr]
Folders are linked (junction on Windows, symlink elsewhere), not copied; nothing under --corpora is modified.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "reference_engine_code"
NAMES = ["breast", "lung", "colorectal", "prostate", "melanoma", "bladder", "kidney", "pancreas", "liver", "stomach", "esophageal",
         "ovarian", "endometrial", "cervical", "thyroid", "brain", "leukemia", "lymphoma", "myeloma"]
CORRECTED = {"lung": "lung_cancer_inclusive_all", "thyroid": "thyroid_cancer_inclusive_all", "kidney": "kidney_cancer_inclusive_all"}
LEXICONS = ["Human Genes.csv", "kegg_compounds_lowercase.csv", "kegg_drug_list_lowercase.csv", "phytochemicals.csv"]


def link(dst: Path, target: Path) -> None:
    if dst.exists() or dst.is_symlink():
        if dst.resolve() == target.resolve():
            return
        raise SystemExit(f"{dst} already exists and points elsewhere")
    dst.parent.mkdir(parents=True, exist_ok=True)
    if sys.platform == "win32":
        r = subprocess.run(["cmd", "/c", "mklink", "/J", str(dst), str(target)], capture_output=True, text=True)
        if r.returncode:
            raise SystemExit(f"could not create junction {dst} -> {target}: {r.stdout}{r.stderr}")
    else:
        dst.symlink_to(target, target_is_directory=True)


def engine_files(root: Path) -> None:
    root.mkdir(parents=True, exist_ok=True)
    for n in LEXICONS:
        shutil.copy2(REF / n, root / n)
    dst = root / "src" / "oncodigger" / "search"
    if not dst.exists():
        shutil.copytree(REF / "oncodigger_search", dst, ignore=shutil.ignore_patterns("__pycache__"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--corpora", type=Path)
    ap.add_argument("--work", type=Path)
    ap.add_argument("--engine-into", type=Path, nargs="+", help="only copy the lexicons and the frozen engine into these (existing) root folders")
    ap.add_argument("--raw", type=Path, help="folder holding <collection>_cancer_inclusive/abstracts.jsonl (optional; see docs/REPRODUCTION.md)")
    ap.add_argument("--corr", action="store_true", help="also build the corrected-collections root (after vocab_corrected_rerun.py)")
    a = ap.parse_args()
    if a.engine_into:
        for r in a.engine_into:
            engine_files(r)
            print("engine files in", r)
        return
    if not (a.corpora and a.work):
        ap.error("--corpora and --work are required")
    a.work.mkdir(parents=True, exist_ok=True)
    frozen = a.work / "s27frozen"
    engine_files(frozen)
    for n in NAMES:
        d = f"{n}_cancer_inclusive_all"
        if not (a.corpora / d / "lexical_index.db").exists():
            raise SystemExit(f"missing {a.corpora / d / 'lexical_index.db'}")
        link(frozen / "data" / "pubmed_processed" / d, a.corpora / d)
    if a.raw:
        link(frozen / "data" / "pubmed_raw", a.raw)
    print("frozen root:", frozen)
    if a.corr:
        corr = a.work / "s27corr"
        engine_files(corr)
        cdir = ROOT / "data" / "derived" / "intogen_rescoring" / "corrected_corpora"
        for n in NAMES:
            d = f"{n}_cancer_inclusive_all"
            src = cdir / d if n in CORRECTED else a.corpora / d
            if not (src / "lexical_index.db").exists():
                raise SystemExit(f"missing {src / 'lexical_index.db'} (run vocab_corrected_rerun.py first)")
            link(corr / "data" / "pubmed_processed" / d, src)
        print("corrected root:", corr)


if __name__ == "__main__":
    main()
