# Environment for the reproduction scripts. Source it from the repository (. scripts/reproduce/env.sh); nothing is hard-coded to one machine.
#
#   ONCODIGGER_CORPORA  (required) folder holding the 19 <collection>_cancer_inclusive_all/lexical_index.db folders of the 2026-09-27 build
#                       (and, if present, pooled_common_cancers_all). SHA-256 of each database: data/provenance/corpus_build_2026-09-27.json
#   REPRO_WORK          (optional) work folder for roots, caches and logs; default <repo>/../oncodigger_repro_work
#   BGFIX_WORKERS       (optional) processes for the exact background count; default 8
#   REPRO_NO_CORPORA=1  (optional) for the corpus-free stages (50, 60), which need no corpus folder
#   REPRO_RESUME=1      (optional) skip steps already recorded as exit=0 in $WORK/logs/progress.txt (continue after a failure)
#   PYTHON              (optional) interpreter; default python
_here="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
REPO="$(cd "$_here/../.." && (pwd -W 2>/dev/null || pwd))"
[ -n "$ONCODIGGER_CORPORA" ] || [ "$REPRO_NO_CORPORA" = 1 ] || { echo "set ONCODIGGER_CORPORA to the folder holding the 19 corpus databases" >&2; return 1 2>/dev/null || exit 1; }
WORK="${REPRO_WORK:-$REPO/../oncodigger_repro_work}"; mkdir -p "$WORK"; WORK="$(cd "$WORK" && (pwd -W 2>/dev/null || pwd))"
PYTHON="${PYTHON:-python}"
export PYTHONHASHSEED=0 PYTHONUTF8=1
export ONCODIGGER_ROOT="$WORK/s27frozen"
export GENE_QUALITY_TEST_ENGINE="$REPO/reference_engine_code/gene_quality_test_engine"
export BGFIX_ENABLE=1
export BGFIX_CACHE="$WORK/bgfix_cache"; mkdir -p "$BGFIX_CACHE"
export BGFIX_WORKERS="${BGFIX_WORKERS:-8}"
export PYTHONPATH="$REPO/exploratory/background_rebuild/hook${PYTHONPATH:+:$PYTHONPATH}"
LOGS="$WORK/logs"; mkdir -p "$LOGS"
cd "$REPO"
# run <name> <script> [args]: runs a script with the hook, keeps its log, records exit code and seconds, stops the stage on failure
run() {
  _n=$1; shift
  if [ "$REPRO_RESUME" = 1 ] && grep -q "^$_n exit=0 " "$LOGS/progress.txt" 2>/dev/null; then echo "$_n skipped (already done)"; return 0; fi
  _t0=$(date +%s)
  BGFIX_STATS="$LOGS/stats_$_n.json" "$PYTHON" -X utf8 "$@" > "$LOGS/$_n.log" 2>&1
  _rc=$?
  echo "$_n exit=$_rc seconds=$(( $(date +%s) - _t0 ))" | tee -a "$LOGS/progress.txt"
  [ $_rc -eq 0 ] || { echo "FAILED: $_n (see $LOGS/$_n.log)" >&2; exit $_rc; }
}
