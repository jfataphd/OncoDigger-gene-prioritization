#!/bin/bash
# Stage 40: post hoc sensitivity analysis without /secondary-indexed records (all 19 collections, two batches). The PMID lists are committed
# (data/derived/metastatic_exclusion/secondary_pmids_*.csv); fetch_secondary_pmids.py (live PubMed query) is not rerun here.
. "$(dirname "$0")/env.sh" || exit 1
NM="$WORK/nometa"; mkdir -p "$NM"
"$PYTHON" -X utf8 scripts/reproduce/assemble_roots.py --engine-into "$NM" || exit 1
B1="liver brain lung kidney melanoma"
B2="breast colorectal prostate bladder pancreas stomach esophageal ovarian endometrial cervical thyroid leukemia lymphoma myeloma"
run m01_build_excluded scripts/analysis/build_metastasis_excluded_corpora.py --src-root "$ONCODIGGER_ROOT" --out-root "$NM" --corpora $B1
run m02_verify_idf scripts/analysis/build_date_truncated_corpora.py --src-root "$NM" --verify --corpora $B1
run m03_rank_excluded scripts/analysis/run_query_planes.py --oncodigger-root "$NM" --allow-nonstandard-build --planes mutation --cancers "Liver Cancer" "Brain Cancer" "Lung Cancer" "Kidney Cancer" "Melanoma" --out-dir "$REPO/data/derived/metastatic_exclusion/ranking"
run n01_build_excluded14 scripts/analysis/build_metastasis_excluded_corpora.py --src-root "$ONCODIGGER_ROOT" --out-root "$NM" --corpora $B2
run n02_verify_idf14 scripts/analysis/build_date_truncated_corpora.py --src-root "$NM" --verify --corpora $B2
run n03_rank_excluded14 scripts/analysis/run_query_planes.py --oncodigger-root "$NM" --allow-nonstandard-build --planes mutation --cancers "Breast Cancer" "Colorectal Cancer" "Prostate Cancer" "Bladder Cancer" "Pancreatic Cancer" "Stomach Cancer" "Esophageal Cancer" "Ovarian Cancer" "Endometrial Cancer" "Cervical Cancer" "Thyroid Cancer" "Leukemia" "Lymphoma" "Myeloma" --out-dir "$REPO/data/derived/metastatic_exclusion/ranking_other14"
echo "STAGE 40 DONE"
