#!/bin/bash
# Stage 45 (run BEFORE stage 50, which reads its outputs): the comparator top-10 lists and their IntOGen scoring.
#   CancerMine (fixed Zenodo record, fetched and checksum-verified), DISEASES 2.0 and PubTator 3.0 (cancer subsets: data/external/diseases and
#   data/external/pubtator; SHA-256 in data/external/PROVENANCE.md; supplied with the Zenodo archive of the release, not committed to GitHub), ChatGPT (the saved response).
# Not rerun here because they query live services whose answers change: pubtator_benchmark_v2_subtype_aggregated.py,
#   pubtator_v1_single_parent_reconstruction.py, cancer_specific_drivers_pubmed_validation.py, fetch_secondary_pmids.py. Their outputs are committed and checksummed.
export REPRO_NO_CORPORA=1
. "$(dirname "$0")/env.sh" || exit 1
unset ONCODIGGER_ROOT BGFIX_ENABLE
"$PYTHON" -X utf8 scripts/reproduce/fetch_external.py || exit 1
for f in data/external/diseases/textmining_full_cancer_subset.tsv data/external/pubtator/relation_cancer_subset.tsv; do
  [ -f "$f" ] || { echo "MISSING $f (supplied with the Zenodo archive of the release; see data/external/PROVENANCE.md)" >&2; exit 1; }
done
"$PYTHON" -X utf8 - <<'PY' || exit 1
import hashlib, re
from pathlib import Path
prov = Path("data/external/PROVENANCE.md").read_text(encoding="utf-8")
for f in ("data/external/diseases/textmining_full_cancer_subset.tsv", "data/external/pubtator/relation_cancer_subset.tsv"):
    h = hashlib.sha256()
    with open(f, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    assert h.hexdigest() in prov, f"{f}: SHA-256 {h.hexdigest()} is not the one recorded in PROVENANCE.md"
    print("checksum OK", f)
PY
run c01_cancermine_benchmark scripts/analysis/cancermine_benchmark.py
run c02_diseases_benchmark scripts/analysis/diseases_benchmark.py
run c03_chatgpt_benchmark scripts/analysis/chatgpt_benchmark.py
run c04_intogen_diseases_pubtator_rescoring scripts/analysis/intogen_diseases_pubtator_rescoring.py
run c05_pool_size_jaccard_unsaturated scripts/analysis/pool_size_jaccard_unsaturated.py
run c06_vocab_policy_v5_table scripts/analysis/vocab_correction_v5_table.py
echo "STAGE 45 DONE"
