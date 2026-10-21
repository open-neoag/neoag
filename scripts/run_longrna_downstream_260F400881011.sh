#!/usr/bin/env bash
set -euo pipefail

# Backward-compatible wrapper for the historical Xiaoliang Site 1 run.
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
WORK=/mnt/zjl-bgi-zzb/peixunban/gl/data/chenxiaoliang_data/work/longrna/260F400881011_20260722

exec "$ROOT/scripts/run_longrna_downstream.sh" \
  --sample-id 260F400881011 \
  --workdir "$WORK" \
  --reference-fasta /root/neo/neodata4git/data/ref/hg38/Homo_sapiens_assembly38.fasta \
  --annotation-gtf /root/neo/neodata4git/data/ref/hg38/gencode.gtf \
  --sqanti-home /root/neo/env_tool/tools/SQANTI3 \
  --sqanti-python /root/neo/env_tool/miniforge3/envs/neoag-sqanti3/bin/python \
  --star-junctions /mnt/zjl-bgi-zzb/peixunban/gl/data/chenxiaoliang_data/work/snaf/FP500004780_L01_203_20260721/star/FP500004780_L01_203.SJ.out.tab \
  --threads "${LONGRNA_DOWNSTREAM_THREADS:-4}" \
  "$@"
