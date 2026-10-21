#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SAMPLE_ID=""; INPUT_DIR=""; WORKDIR=""; REFERENCE_FASTA=""; REFERENCE_INDEX=""
ANNOTATION_GTF=""; HLA_FILE=""; SHORTREAD_COMPARISONS=""; NORMAL_JUNCTIONS=""
ILLUMINA_BAM=""; STAR_JUNCTIONS=""; ISOQUANT_BIN="isoquant"; SQANTI_HOME=""
SQANTI_PYTHON="python3"; THREADS=8; PEPTIDE_LENGTHS="8,9,10,11"
RETENTION="standard"; RESUME=0; DRY_RUN=0; PATCH_SQANTI3=0

usage() {
  cat <<'EOF'
Usage: run_longrna_splice_full.sh [options]

Required:
  --sample-id ID                 Sample identifier propagated to every output
  --input-dir DIR                Directory containing long-read FASTQ files
  --workdir DIR                  Canonical output directory
  --reference-fasta FILE         GRCh38 reference FASTA
  --annotation-gtf FILE          Matching gene annotation GTF
  --hla FILE                     One HLA allele per line
  --shortread-comparisons FILE   TSV containing junction_key values
  --sqanti-home DIR              SQANTI3 installation containing sqanti3_qc.py

Options:
  --reference-index FILE         Existing minimap2 index for IsoQuant
  --illumina-bam FILE            Optional matched short-read RNA BAM for IsoQuant
  --star-junctions FILE          Optional STAR SJ.out.tab for SQANTI3 coverage
  --normal-junctions FILE        Optional normal-cohort junction_key TSV
  --isoquant-bin FILE            IsoQuant executable (default: isoquant)
  --sqanti-python FILE           SQANTI3 Python interpreter (default: python3)
  --threads N                    Worker count (default: 8)
  --peptide-lengths CSV          MHC-I peptide lengths (default: 8,9,10,11)
  --retention MODE               standard or full (default: standard)
  --patch-sqanti3                Apply the audited SQANTI3 compatibility patch
  --resume                       Reuse stages with validated completion markers
  --dry-run                      Validate inputs and print the execution plan
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-id) SAMPLE_ID=$2; shift 2 ;;
    --input-dir) INPUT_DIR=$2; shift 2 ;;
    --workdir) WORKDIR=$2; shift 2 ;;
    --reference-fasta) REFERENCE_FASTA=$2; shift 2 ;;
    --reference-index) REFERENCE_INDEX=$2; shift 2 ;;
    --annotation-gtf) ANNOTATION_GTF=$2; shift 2 ;;
    --hla) HLA_FILE=$2; shift 2 ;;
    --shortread-comparisons) SHORTREAD_COMPARISONS=$2; shift 2 ;;
    --normal-junctions) NORMAL_JUNCTIONS=$2; shift 2 ;;
    --illumina-bam) ILLUMINA_BAM=$2; shift 2 ;;
    --star-junctions) STAR_JUNCTIONS=$2; shift 2 ;;
    --isoquant-bin) ISOQUANT_BIN=$2; shift 2 ;;
    --sqanti-home) SQANTI_HOME=$2; shift 2 ;;
    --sqanti-python) SQANTI_PYTHON=$2; shift 2 ;;
    --threads) THREADS=$2; shift 2 ;;
    --peptide-lengths) PEPTIDE_LENGTHS=$2; shift 2 ;;
    --retention) RETENTION=$2; shift 2 ;;
    --patch-sqanti3) PATCH_SQANTI3=1; shift ;;
    --resume) RESUME=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$SAMPLE_ID" && -n "$INPUT_DIR" && -n "$WORKDIR" && -n "$REFERENCE_FASTA" && -n "$ANNOTATION_GTF" ]] || { usage >&2; exit 2; }
[[ -n "$HLA_FILE" && -n "$SHORTREAD_COMPARISONS" && -n "$SQANTI_HOME" ]] || { usage >&2; exit 2; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --threads must be a positive integer" >&2; exit 2; }
[[ "$RETENTION" == standard || "$RETENTION" == full ]] || { echo "ERROR: --retention must be standard or full" >&2; exit 2; }

required=("$INPUT_DIR" "$REFERENCE_FASTA" "$ANNOTATION_GTF" "$HLA_FILE" "$SHORTREAD_COMPARISONS" "$SQANTI_HOME/sqanti3_qc.py")
for path in "${required[@]}"; do [[ -e "$path" ]] || { echo "ERROR: required path not found: $path" >&2; exit 1; }; done
for path in "$REFERENCE_INDEX" "$ILLUMINA_BAM" "$STAR_JUNCTIONS" "$NORMAL_JUNCTIONS"; do
  [[ -z "$path" || -s "$path" ]] || { echo "ERROR: optional input is missing or empty: $path" >&2; exit 1; }
done
mapfile -t FASTQS < <(find "$INPUT_DIR" -maxdepth 1 -type f \( -name '*.fastq' -o -name '*.fastq.gz' -o -name '*.fq' -o -name '*.fq.gz' \) -print | sort)
(( ${#FASTQS[@]} > 0 )) || { echo "ERROR: no FASTQ files found in $INPUT_DIR" >&2; exit 1; }

isoquant_args=(
  --reference "$REFERENCE_FASTA"
  --genedb "$ANNOTATION_GTF"
  --complete_genedb
  --fastq "${FASTQS[@]}"
  --data_type nanopore
  --stranded none
  --threads "$THREADS"
  --check_canonical
  --sqanti_output
  --count_exons
  --read_group none
  --large_output read_assignments corrected_bed read2transcripts
  --output "$WORKDIR/isoquant"
  --prefix "$SAMPLE_ID"
)
[[ -n "$REFERENCE_INDEX" ]] && isoquant_args+=(--index "$REFERENCE_INDEX")
[[ -n "$ILLUMINA_BAM" ]] && isoquant_args+=(--illumina_bam "$ILLUMINA_BAM")
[[ "$RETENTION" == full ]] && isoquant_args+=(--keep_tmp)

if (( DRY_RUN )); then
  printf 'sample_id=%s\nfastq_files=%s\nworkdir=%s\nretention=%s\n' "$SAMPLE_ID" "${#FASTQS[@]}" "$WORKDIR" "$RETENTION"
  printf 'isoquant_command='; printf '%q ' "$ISOQUANT_BIN" "${isoquant_args[@]}"; printf '\n'
  printf 'next_stages=%s\n' 'SQANTI3+TD2,candidate_filter,normal_background,MHC_prediction,evidence_merge,verification'
  [[ -n "$NORMAL_JUNCTIONS" ]] || printf 'normal_background=UNASSESSED\n'
  exit 0
fi

if [[ -e "$WORKDIR/longrna_splice_status.tsv" && $RESUME -eq 0 ]]; then
  echo "ERROR: workdir already contains a run; pass --resume or choose a new --workdir" >&2
  exit 1
fi
if [[ ! -e "$WORKDIR/longrna_splice_status.tsv" ]]; then
  "$ROOT/scripts/run_longrna_splice.sh" --sample-id "$SAMPLE_ID" --input-dir "$INPUT_DIR" --workdir "$WORKDIR"
fi
STATUS="$WORKDIR/longrna_splice_status.tsv"
record() { printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$(date -Iseconds)" "$3" >> "$STATUS"; }
mkdir -p "$WORKDIR"/{logs,translation,splice_candidates,presentation,integration,provenance}
COMMAND_LOG="$WORKDIR/provenance/commands.sh"
if [[ ! -e "$COMMAND_LOG" ]]; then
  printf '#!/usr/bin/env bash\nset -euo pipefail\n' > "$COMMAND_LOG"
  chmod +x "$COMMAND_LOG"
fi
log_command() { printf '%q ' "$@" >> "$COMMAND_LOG"; printf '\n' >> "$COMMAND_LOG"; }

provenance_inputs=(
  "$SQANTI_PYTHON" "$ROOT/scripts/capture_longrna_provenance.py"
  --sample-id "$SAMPLE_ID"
  --retention "$RETENTION"
  --scope input
  --output "$WORKDIR/provenance/input_manifest.json"
  --git-root "$ROOT"
  --path "$INPUT_DIR"
  --path "$REFERENCE_FASTA"
  --path "$ANNOTATION_GTF"
  --path "$HLA_FILE"
  --path "$SHORTREAD_COMPARISONS"
  --tool "isoquant=$ISOQUANT_BIN"
  --tool "sqanti_python=$SQANTI_PYTHON"
)
for path in "$REFERENCE_INDEX" "$ILLUMINA_BAM" "$STAR_JUNCTIONS" "$NORMAL_JUNCTIONS"; do
  [[ -n "$path" ]] && provenance_inputs+=(--path "$path")
done
log_command "${provenance_inputs[@]}"
"${provenance_inputs[@]}"
record provenance_capture PASS "retention=$RETENTION input_manifest=$WORKDIR/provenance/input_manifest.json"

RUN_COMPLETE=0
on_exit() {
  code=$?
  if (( code != 0 && RUN_COMPLETE == 0 )); then
    record workflow FAILED "exit_code=$code"
    failure_args=(
      "$SQANTI_PYTHON" "$ROOT/scripts/capture_longrna_provenance.py"
      --sample-id "$SAMPLE_ID" --retention "$RETENTION" --scope failure
      --output "$WORKDIR/provenance/failure_manifest.json"
      --path "$WORKDIR/logs" --path "$STATUS" --path "$COMMAND_LOG"
    )
    [[ "$RETENTION" == full ]] && failure_args+=(--path "$WORKDIR")
    "${failure_args[@]}" >/dev/null 2>&1 || true
  fi
  exit "$code"
}
trap on_exit EXIT

if [[ -e "$WORKDIR/.isoquant_complete" ]]; then
  mapfile -t gtfs < <(find "$WORKDIR/isoquant" -type f -name '*.transcript_models.gtf' -size +0c -print)
  (( ${#gtfs[@]} == 1 )) || {
    record isoquant_transcript_reconstruction FAILED "completion marker exists but transcript_model_gtf_count=${#gtfs[@]}"
    exit 1
  }
  record isoquant_transcript_reconstruction REUSED "validated_gtf=${gtfs[0]}"
else
  record isoquant_transcript_reconstruction RUNNING "fastq_files=${#FASTQS[@]}"
  log_command "$ISOQUANT_BIN" "${isoquant_args[@]}"
  if "$ISOQUANT_BIN" "${isoquant_args[@]}" >> "$WORKDIR/logs/isoquant.log" 2>&1; then
    mapfile -t gtfs < <(find "$WORKDIR/isoquant" -type f -name '*.transcript_models.gtf' -size +0c -print)
    (( ${#gtfs[@]} == 1 )) || { record isoquant_transcript_reconstruction FAILED "transcript_model_gtf_count=${#gtfs[@]}"; exit 1; }
    touch "$WORKDIR/.isoquant_complete"
    record isoquant_transcript_reconstruction PASS "gtf=${gtfs[0]}"
  else
    code=$?; record isoquant_transcript_reconstruction FAILED "exit_code=$code log=$WORKDIR/logs/isoquant.log"; exit "$code"
  fi
fi

downstream=(
  "$ROOT/scripts/run_longrna_downstream.sh"
  --sample-id "$SAMPLE_ID"
  --workdir "$WORKDIR"
  --reference-fasta "$REFERENCE_FASTA"
  --annotation-gtf "$ANNOTATION_GTF"
  --sqanti-home "$SQANTI_HOME"
  --sqanti-python "$SQANTI_PYTHON"
  --threads "$THREADS"
)
[[ -n "$STAR_JUNCTIONS" ]] && downstream+=(--star-junctions "$STAR_JUNCTIONS")
(( PATCH_SQANTI3 )) && downstream+=(--patch-sqanti3)
log_command "${downstream[@]}"
"${downstream[@]}"

protein=$(<"$WORKDIR/sqanti3/protein_output.path")
ln -sfn "$protein" "$WORKDIR/translation/${SAMPLE_ID}.pep"
candidate_args=(
  --sample-id "$SAMPLE_ID"
  --sqanti-classification "$WORKDIR/sqanti3/${SAMPLE_ID}_classification.txt"
  --sqanti-junctions "$WORKDIR/sqanti3/${SAMPLE_ID}_junctions.txt"
  --corrected-gtf "$WORKDIR/sqanti3/${SAMPLE_ID}_corrected.gtf"
  --td2-peptides "$protein"
  --shortread-comparisons "$SHORTREAD_COMPARISONS"
  --hla "$HLA_FILE"
  --outdir "$WORKDIR/splice_candidates"
  --peptide-lengths "$PEPTIDE_LENGTHS"
)
[[ -n "$NORMAL_JUNCTIONS" ]] && candidate_args+=(--normal-junctions "$NORMAL_JUNCTIONS")
log_command "$SQANTI_PYTHON" "$ROOT/scripts/build_longrna_unique_splice_peptides.py" "${candidate_args[@]}"
"$SQANTI_PYTHON" "$ROOT/scripts/build_longrna_unique_splice_peptides.py" "${candidate_args[@]}"
record junction_peptide_generation PASS "outdir=$WORKDIR/splice_candidates"

predict_command=(bash "$ROOT/scripts/run_peptide_predict.sh" \
  --input "$WORKDIR/splice_candidates/longrna_unique_prediction_input.tsv" \
  --outdir "$WORKDIR/presentation" \
  --sample-id "$SAMPLE_ID")
log_command "${predict_command[@]}"
"${predict_command[@]}"
record mhc_presentation_prediction PASS "outdir=$WORKDIR/presentation"

mapfile -t presentation_files < <(find "$WORKDIR/presentation" -type f -name 'presentation_evidence.tsv' -size +0c -print)
(( ${#presentation_files[@]} == 1 )) || {
  record comprehensive_evidence FAILED "presentation_evidence_count=${#presentation_files[@]}"
  exit 1
}
merge_command=("$SQANTI_PYTHON" "$ROOT/scripts/merge_longrna_comprehensive_evidence.py" \
  --provenance "$WORKDIR/splice_candidates/longrna_unique_junction_peptides.tsv" \
  --presentation "${presentation_files[0]}" \
  --output "$WORKDIR/comprehensive_evidence.tsv" \
  --summary "$WORKDIR/integration/comprehensive_evidence.summary.json")
log_command "${merge_command[@]}"
"${merge_command[@]}"
record comprehensive_evidence PASS "output=$WORKDIR/comprehensive_evidence.tsv"

log_command bash "$ROOT/scripts/verify_longrna_splice.sh" "$WORKDIR"
bash "$ROOT/scripts/verify_longrna_splice.sh" "$WORKDIR" | tee "$WORKDIR/integration/verification.txt"
record workflow PASS "retention=$RETENTION"

output_manifest=(
  "$SQANTI_PYTHON" "$ROOT/scripts/capture_longrna_provenance.py"
  --sample-id "$SAMPLE_ID" --retention "$RETENTION" --scope output
  --output "$WORKDIR/provenance/output_manifest.json" --path "$WORKDIR"
)
log_command "${output_manifest[@]}"
"${output_manifest[@]}"
RUN_COMPLETE=1
trap - EXIT
