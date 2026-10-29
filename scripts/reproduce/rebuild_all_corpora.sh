#!/bin/bash
# Download the 19 collections of the 2026-09-27 build from NCBI by PMID and rebuild their lexical_index.db files.
# Not part of run_all.sh: it downloads about 1.6 million PubMed records (hours; much faster with an NCBI API key) and is meant for readers who do not hold
# the author's databases. Resumes where it stopped (raw downloads are cached). After it finishes, point ONCODIGGER_CORPORA at <out>/pubmed_processed
# and run the stages. A rebuild can differ slightly from the databases used in the paper (PubMed records are revised); see docs/REPRODUCTION.md.
#   usage: bash scripts/reproduce/rebuild_all_corpora.sh <out-folder> [NCBI_API_KEY]
set -e
OUT="${1:?usage: rebuild_all_corpora.sh <out-folder> [NCBI_API_KEY]}"; KEY="${2:-}"
cd "$(dirname "$0")/../.."
for c in $(tail -n +2 data/provenance/pmid_manifests/MANIFEST.csv | cut -d, -f1); do
  d="$OUT/pubmed_processed/${c}_cancer_inclusive_all"; mkdir -p "$d"
  [ -f "$d/lexical_index.db" ] && { echo "== $c already built"; continue; }
  echo "== $c"
  python -X utf8 scripts/reproduce/rebuild_corpus_from_pmids.py "$c" --out "$OUT/work/$c" ${KEY:+--api-key "$KEY"}
  mv "$OUT/work/$c/${c}_lexical_index.db" "$d/lexical_index.db"
done
echo "Done. Set ONCODIGGER_CORPORA=$OUT/pubmed_processed"
