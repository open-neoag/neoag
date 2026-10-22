#!/usr/bin/env bash
set -euo pipefail
WORKDIR=${1:?Usage: verify_longrna_splice.sh WORKDIR}
STATUS="$WORKDIR/longrna_splice_status.tsv"
fail=0
check() {
  if [[ -s "$1" ]]; then printf 'PASS\t%s\n' "$1"; else printf 'MISSING\t%s\n' "$1"; fail=1; fi
}
check "$STATUS"
check "$WORKDIR/README.workflow.txt"
if [[ -s "$WORKDIR/isoquant_transcript_model.path" ]]; then
  isoform_gtf=$(<"$WORKDIR/isoquant_transcript_model.path")
  if [[ -s "$isoform_gtf" ]]; then
    printf 'PASS\t%s -> %s\n' "$WORKDIR/isoquant_transcript_model.path" "$isoform_gtf"
  else
    printf 'MISSING\t%s -> %s\n' "$WORKDIR/isoquant_transcript_model.path" "$isoform_gtf"
    fail=1
  fi
else
  found=0
  for f in "$WORKDIR"/isoquant/*/*.transcript_models.gtf; do
    if [[ -e "$f" ]]; then check "$f"; found=1; fi
  done
  if (( found == 0 )); then printf 'MISSING\t%s\n' "$WORKDIR/isoquant/*/*.transcript_models.gtf"; fail=1; fi
fi
for pattern in \
  "$WORKDIR/sqanti3/*_classification.txt" \
  "$WORKDIR/translation/*.pep" \
  "$WORKDIR/comprehensive_evidence.tsv"; do
  found=0
  for f in $pattern; do
    if [[ -e "$f" ]]; then check "$f"; found=1; fi
  done
  if (( found == 0 )); then printf 'MISSING\t%s\n' "$pattern"; fail=1; fi
done
if awk -F '\t' 'NR > 1 {latest[$1]=$2} END {for (stage in latest) if (latest[stage] == "FAILED") exit 0; exit 1}' "$STATUS"; then
  echo 'FAILED\tlatest stage status contains a failure'
  fail=1
fi
if (( fail == 0 )); then echo workflow_status=PASS; else echo workflow_status=INCOMPLETE; exit 1; fi
