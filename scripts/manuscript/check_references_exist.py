"""
check_references_exist.py -- every file path and every supplementary label the manuscript cites must exist.

* Backticked repository paths (scripts/..., data/..., docs/..., figures/..., manuscript/...) must exist in the repository.
  Backticked bare file names ending .py or .csv must match a file somewhere in the repository.
* Supplementary Tables S<n> and repository tables RT<n> cited anywhere must be defined (a "**Supplementary Table S<n>." heading, or an
  RT<n> file in data/supplementary_tables/); Supplementary Figures and Supplementary Results/Methods S<n> likewise.
Exit 1 and a list if anything is missing. Run: python -X utf8 scripts/manuscript/check_references_exist.py [manuscript.md]
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MD = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "manuscript" / "OncoDigger_manuscript1.1.md"
TEXT = MD.read_text(encoding="utf-8")
SKIP = {".git", "node_modules", "__pycache__"}
files = [p for p in ROOT.rglob("*") if p.is_file() and not (SKIP & set(p.relative_to(ROOT).parts))]
names = {p.name for p in files}
problems: list[str] = []

for m in sorted(set(re.findall(r"`([^`\s]+)`", TEXT))):
    path = m.rstrip(".,;:)")
    if re.match(r"^(scripts|data|docs|figures|manuscript|reference_engine_code|exploratory)/", path):
        if "*" in path or "<" in path:
            continue
        if not (ROOT / path.rstrip("/")).exists():
            problems.append(f"missing path: {path}")
    elif re.match(r"^[\w\-.]+\.(py|csv)$", path) and path not in names:
        problems.append(f"missing file name: {path}")

defined_t = set(re.findall(r"\*\*Supplementary Table (S\d+)\b", TEXT))
cited_t = set(re.findall(r"Supplementary Tables? (S\d+)\b", TEXT)) | set(re.findall(r"\b(S\d+)\b(?= and S\d+| to S\d+)", TEXT))
for t in sorted(cited_t - defined_t, key=lambda x: int(x[1:])):
    problems.append(f"Supplementary Table {t} cited but not defined")

rt_files = {m.group(1) for p in (ROOT / "data" / "supplementary_tables").glob("RT*") if (m := re.match(r"(RT\d+)[AB]?_", p.name))}
for t in sorted(set(re.findall(r"\b(RT\d+)\b", TEXT)) - rt_files, key=lambda x: int(x[2:])):
    problems.append(f"{t} cited but no file in data/supplementary_tables/")

for kind in ("Supplementary Figure", "Supplementary Results", "Supplementary Methods"):
    d = set(re.findall(r"^#+ " + kind + r" (S\d+)", TEXT, re.M))
    c = set(re.findall(kind + r" (S\d+)", TEXT))
    for t in sorted(c - d, key=lambda x: int(x[1:])):
        problems.append(f"{kind} {t} cited but not defined")

print(f"{len(files)} repository files; {len(defined_t)} supplementary tables defined; {len(rt_files)} repository tables")
if problems:
    print(f"{len(problems)} problem(s):")
    for x in problems:
        print("  ", x)
    sys.exit(1)
print("all cited paths and labels exist")
