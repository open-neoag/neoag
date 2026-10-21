# Long-read RNA splice workflow

`scripts/run_longrna_splice_full.sh` is the supported end-to-end entry point for a new sample. It runs IsoQuant, SQANTI3 with TD2 ORF prediction, long-read-specific junction filtering, peptide generation, presentation prediction, evidence merging, and output verification.

## Required inputs

- Long-read FASTQ files from one sample
- A genome FASTA and matching annotation GTF
- A sample HLA allele list
- A TSV with a `junction_key` column for short-read comparison
- A SQANTI3 installation and compatible Python environment

A normal-cohort TSV with a `junction_key` column is optional. If it is omitted, normal-junction evidence is recorded as `UNASSESSED`; absence of that input is never interpreted as a negative normal-background result.

## Example

```bash
bash scripts/run_longrna_splice_full.sh \
  --sample-id CASE001 \
  --input-dir /data/CASE001/longrna_fastq \
  --workdir /work/CASE001/longrna \
  --reference-fasta /refs/GRCh38.fa \
  --reference-index /refs/GRCh38.mmi \
  --annotation-gtf /refs/gencode.gtf \
  --hla /data/CASE001/hla.txt \
  --shortread-comparisons /data/CASE001/shortread_junctions.tsv \
  --normal-junctions /refs/normal_junctions.tsv \
  --isoquant-bin /envs/longrna/bin/isoquant \
  --sqanti-home /tools/SQANTI3 \
  --sqanti-python /envs/sqanti3/bin/python \
  --retention full \
  --threads 16
```

Use `--dry-run` to validate paths and inspect the IsoQuant command. Use `--resume` only with the same sample, references, and input definitions; completed stages are reused through explicit marker files.

## Retention and provenance

`--retention standard` is the default. It preserves the stable stage outputs, logs, command history, runtime/package inventory, and file metadata without hashing every large file.

`--retention full` additionally passes `--keep_tmp` to IsoQuant and computes SHA-256 for every inventoried input and output file. It retains the workflow work directory on failure and writes a failure manifest. Full retention can consume substantial storage and the initial/final hashing passes can take significant time for large FASTQ, BAM, and temporary files.

Both modes write reproducibility records under `provenance/`:

- `commands.sh`: shell-escaped stage commands in execution order
- `input_manifest.json`: inputs, Git state, runtime, tools, Python packages, and Conda packages
- `output_manifest.json`: completed output inventory
- `failure_manifest.json`: logs and failure-state inventory when a stage exits unsuccessfully

Raw inputs and shared references remain in their original locations. Full retention records their resolved paths and SHA-256 values; it does not duplicate them into the work directory.

`scripts/run_longrna_downstream.sh` is the parameterized SQANTI3/TD2 stage entry point. The historical `run_longrna_downstream_260F400881011.sh` remains as a compatibility wrapper for Xiaoliang Site 1.

The compatibility patch and recovery utilities are not routine stages:

- `fix_sqanti3_isoform_hits_api.py` patches the audited historical SQANTI3 API conflict.
- `resume_sqanti3_outputs.py` reconstructs outputs from the corresponding historical partial run.
