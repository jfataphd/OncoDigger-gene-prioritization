"""
check_citations.py — cross-check in-text citations against a manuscript's own References list.

Reports (1) in-text citations with no matching reference (first-author surname + year, including a/b suffixes)
and (2) references that are never cited in the text.

Defaults to OncoDigger_manuscript1.1.md, the actual submission (this script previously hardcoded
OncoDigger_reference.md -- the non-submitted reference document -- which meant the real submission's
citations were never checked by this tool; fixed 2026-10-02). Pass a path to check a different file.

The reference-line regex also previously required a literal ", YYYY." immediately after the author
list (comma then year then period), but this manuscript's AMA/Vancouver style has a journal-title
period between the authors and the year, and the year itself is followed by a semicolon, not a
period (e.g. "...Nucleic Acids Res. 2017;45(D1):D777-D783."). That comma requirement never matched
a single reference line in either document -- this tool has likely always reported "0 references"
regardless of target file. Fixed 2026-10-02 to just require a 4-digit year followed by "." or ";"
anywhere after the author list, not specifically after a comma.

Run:  python -X utf8 scripts/manuscript/check_citations.py [path/to/manuscript.md]
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

MD = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[2] / "manuscript" / "OncoDigger_manuscript1.1.md"
NAME = r"[A-Z][A-Za-zÀ-ÿ'\-]+(?: [A-Z][A-Za-zÀ-ÿ'\-]+)?"
CITE = re.compile(rf"({NAME})(?: et al\.| and {NAME})?,? ((?:19|20)\d\d[a-z]?)")


def main() -> None:
    text = MD.read_text(encoding="utf-8")
    body, rest = text.split("## References", 1)
    refs_block = rest.split("\n## ", 1)[0]
    tail = rest.split("\n## ", 1)[1] if "\n## " in rest else ""
    refs = {}
    for line in refs_block.strip().splitlines():
        line = line.strip()
        m = re.match(r"([A-Z][A-Za-zÀ-ÿ'\-]+).*?((?:19|20)\d\d[a-z]?)[.;]", line)
        if m:
            refs[(m.group(1), m.group(2))] = line[:90]
    cited = set()
    for chunk in re.findall(r"\(([^()]*(?:19|20)\d\d[a-z]?[^()]*)\)", body + tail):
        for m in CITE.finditer(chunk):
            if m.group(1).split()[-1] in {"May", "NDA", "BLA", "January", "June"}:
                continue
            cited.add((m.group(1), m.group(2)))
    surname_year = {(k[0].split()[-1], k[1]) for k in refs}
    missing = sorted(c for c in cited if (c[0].split()[-1], c[1]) not in surname_year)
    unused = sorted(k for k in refs if not any(k[0].split()[-1] == c[0].split()[-1] and k[1] == c[1] for c in cited))
    print(f"{len(refs)} references, {len(cited)} distinct in-text citations")
    print("Cited but not in References:", missing or "none")
    print("In References but never cited:")
    for k in unused:
        print("  ", refs[k])
    if missing or unused:
        # Warning only: the flags are known matcher limits (GeneReviews chapters cited by their book date, "Research Network" for
        # The Cancer Genome Atlas Research Network, month names read as surnames) and were checked by hand; a reviewer reads this list.
        print(f"WARNING: {len(missing)} in-text citation(s) not matched to a reference and {len(unused)} reference(s) not matched to a citation (see above).")
    sys.exit(0)


if __name__ == "__main__":
    main()
