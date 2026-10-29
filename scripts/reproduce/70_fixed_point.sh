#!/bin/bash
# Stage 70: acceptance check. Reruns the downstream stage (50) and verifies that nothing it writes changes: if a result depended on the order
# in which stages ran, the second pass would differ from the first. Then lists every tracked file that differs from the commit.
export REPRO_NO_CORPORA=1
. "$(dirname "$0")/env.sh" || exit 1
# The snapshot hashes every file stage 50 can write (data/ except the shipped external inputs, and manuscript/); it does not need git, so it also works on the unpacked archive.
snap() { find data manuscript -type f -not -path 'data/external/*' -not -name '*.pyc' | sort | while read -r f; do sha256sum "$f"; done | sha256sum; }
[ "$(find data manuscript -type f -not -path 'data/external/*' | wc -l)" -gt 100 ] || { echo "stage 70: expected files under data/ and manuscript/ are missing; run it from the repository root" >&2; exit 1; }
before=$(snap)
bash scripts/reproduce/50_downstream.sh > "$LOGS/stage_50_second_pass.out" 2>&1 || { echo "second pass of stage 50 FAILED"; exit 1; }
after=$(snap)
if [ "$before" = "$after" ]; then echo "FIXED POINT: a second pass of stage 50 changed nothing"; else echo "NOT A FIXED POINT: the second pass changed files; stage order matters" >&2; exit 1; fi
if git rev-parse --is-inside-work-tree > /dev/null 2>&1; then
  echo "Files that differ from the commit:"; git status --short | grep -v '^??' || echo "(none)"
else
  echo "Not a git checkout, so there is no commit to compare with; stage 60 verifies every recorded output against data/provenance/outputs_record.csv."
fi
