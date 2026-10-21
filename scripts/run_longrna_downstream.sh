#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SAMPLE_ID=""
WORKDIR=""
REFERENCE_FASTA=""
ANNOTATION_GTF=""
SQANTI_HOME=""
SQANTI_PYTHON="python3"
ISOFORM_GTF=""
STAR_JUNCTIONS=""
THREADS=8
PATCH_SQANTI3=0
DRY_RUN=0

usage() {
  cat <<'EOF'
Usage: run_longrna_downstream.sh [options]

Required:
  --sample-id ID
  --workdir DIR
  --reference-fasta FILE
  --annotation-gtf FILE
  --sqanti-home DIR

Options:
  --sqanti-python FILE      Python interpreter for SQANTI3 (default: python3)
  --isoform-gtf FILE        IsoQuant transcript-model GTF; auto-detected otherwise
  --star-junctions FILE     Optional STAR SJ.out.tab short-read support
  --threads N               SQANTI3 worker count (default: 8)
  --patch-sqanti3           Apply the audited isoform-hits compatibility patch
  --dry-run                 Print the resolved command without executing it
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-id) SAMPLE_ID=$2; shift 2 ;;
    --workdir) WORKDIR=$2; shift 2 ;;
    --reference-fasta) REFERENCE_FASTA=$2; shift 2 ;;
    --annotation-gtf) ANNOTATION_GTF=$2; shift 2 ;;
    --sqanti-home) SQANTI_HOME=$2; shift 2 ;;
    --sqanti-python) SQANTI_PYTHON=$2; shift 2 ;;
    --isoform-gtf) ISOFORM_GTF=$2; shift 2 ;;
    --star-junctions) STAR_JUNCTIONS=$2; shift 2 ;;
    --threads) THREADS=$2; shift 2 ;;
    --patch-sqanti3) PATCH_SQANTI3=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$SAMPLE_ID" && -n "$WORKDIR" && -n "$REFERENCE_FASTA" && -n "$ANNOTATION_GTF" && -n "$SQANTI_HOME" ]] || {
  usage >&2
  exit 2
}
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --threads must be a positive integer" >&2; exit 2; }

for path in "$REFERENCE_FASTA" "$ANNOTATION_GTF" "$SQANTI_HOME/sqanti3_qc.py"; do
  [[ -s "$path" ]] || { echo "ERROR: required input is missing or empty: $path" >&2; exit 1; }
done
if [[ -n "$STAR_JUNCTIONS" && ! -s "$STAR_JUNCTIONS" ]]; then
  echo "ERROR: --star-junctions is missing or empty: $STAR_JUNCTIONS" >&2
  exit 1
fi

mkdir -p "$WORKDIR"/{logs,sqanti3}
STATUS="$WORKDIR/longrna_splice_status.tsv"
[[ -e "$STATUS" ]] || printf 'stage\tstatus\tupdated_at\tdetail\n' > "$STATUS"
record() { printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$(date -Iseconds)" "$3" >> "$STATUS"; }

validate_outputs() {
  classification="$WORKDIR/sqanti3/${SAMPLE_ID}_classification.txt"
  junctions="$WORKDIR/sqanti3/${SAMPLE_ID}_junctions.txt"
  corrected_gtf="$WORKDIR/sqanti3/${SAMPLE_ID}_corrected.gtf"
  corrected_fasta="$WORKDIR/sqanti3/${SAMPLE_ID}_corrected.fasta"
  for path in "$classification" "$junctions" "$corrected_gtf" "$corrected_fasta"; do
    [[ -s "$path" ]] || { echo "ERROR: SQANTI3 output is missing or empty: $path" >&2; return 1; }
  done
  protein_candidates=(
    "$WORKDIR/sqanti3/${SAMPLE_ID}_corrected.faa"
    "$WORKDIR/sqanti3/TD2/${SAMPLE_ID}_corrected.fasta.TD2.pep"
  )
  protein=""
  for path in "${protein_candidates[@]}"; do
    if [[ -s "$path" ]]; then protein=$path; break; fi
  done
  [[ -n "$protein" ]] || { echo "ERROR: SQANTI3/TD2 protein output is missing" >&2; return 1; }
  printf '%s\n' "$protein" > "$WORKDIR/sqanti3/protein_output.path"
}

if [[ -z "$ISOFORM_GTF" ]]; then
  mapfile -t candidates < <(find "$WORKDIR/isoquant" -type f -name '*.transcript_models.gtf' -size +0c -print | sort)
  (( ${#candidates[@]} == 1 )) || {
    echo "ERROR: expected one IsoQuant transcript-model GTF, found ${#candidates[@]}; pass --isoform-gtf" >&2
    exit 1
  }
  ISOFORM_GTF=${candidates[0]}
fi
[[ -s "$ISOFORM_GTF" ]] || { echo "ERROR: IsoQuant GTF is missing or empty: $ISOFORM_GTF" >&2; exit 1; }

if (( PATCH_SQANTI3 )); then
  "$SQANTI_PYTHON" "$ROOT/scripts/fix_sqanti3_isoform_hits_api.py" --sqanti-home "$SQANTI_HOME"
fi

args=(
  --isoforms "$ISOFORM_GTF"
  --refGTF "$ANNOTATION_GTF"
  --refFasta "$REFERENCE_FASTA"
  --include_ORF
  --isoform_hits
  --report html
  --cpus "$THREADS"
  --chunks 1
  --output "$SAMPLE_ID"
  --dir "$WORKDIR/sqanti3"
)
[[ -n "$STAR_JUNCTIONS" ]] && args+=(--coverage "$STAR_JUNCTIONS")

if (( DRY_RUN )); then
  printf 'isoform_gtf=%s\n' "$ISOFORM_GTF"
  printf 'command='; printf '%q ' "$SQANTI_PYTHON" "$SQANTI_HOME/sqanti3_qc.py" "${args[@]}"; printf '\n'
  exit 0
fi

if [[ -e "$WORKDIR/.sqanti3_complete" ]]; then
  if ! validate_outputs; then
    record sqanti3_structure_orf_annotation FAILED "completion marker exists but outputs failed validation"
    exit 1
  fi
  record sqanti3_structure_orf_annotation REUSED "validated protein=$protein"
  exit 0
fi

record sqanti3_structure_orf_annotation RUNNING "isoforms=$ISOFORM_GTF"
if "$SQANTI_PYTHON" "$SQANTI_HOME/sqanti3_qc.py" "${args[@]}" >> "$WORKDIR/logs/sqanti3.log" 2>&1; then
  if ! validate_outputs; then
    record sqanti3_structure_orf_annotation FAILED "SQANTI3 completed but outputs failed validation"
    exit 1
  fi
  touch "$WORKDIR/.sqanti3_complete"
  record sqanti3_structure_orf_annotation PASS "protein=$protein"
else
  code=$?
  record sqanti3_structure_orf_annotation FAILED "exit_code=$code log=$WORKDIR/logs/sqanti3.log"
  exit "$code"
fi
