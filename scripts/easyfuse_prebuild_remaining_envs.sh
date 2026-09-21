#!/usr/bin/env bash
# Sequential prebuild for remaining EasyFuse conda envs (mamba -y; do not rely on Nextflow).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck source=/dev/null
source "${ROOT}/conf/tools.env.sh"

RUNTIME_DIR="${EASYFUSE_RUNTIME_DIR:-${ROOT}/work}"
CONDA_CACHE="${EASYFUSE_NXF_CONDA_CACHEDIR:-${RUNTIME_DIR}/.nextflow_conda}"
LOG="${EASYFUSE_PREBUILD_LOG:-${RUNTIME_DIR}/easyfuse_conda_prebuild.log}"
mkdir -p "${CONDA_CACHE}" "$(dirname "${LOG}")"

export CONDA_ALWAYS_YES=1

wait_for_mamba_free() {
  while pgrep -f 'mamba env create' >/dev/null 2>&1; do
    sleep 15
  done
  rm -f "${NEOAG_CONDA_BASE}/pkgs/pkgs.lock" 2>/dev/null || true
}

prebuild_conda_env() {
  local env_id="$1"
  local yml="$2"
  local check_bin="$3"
  local prefix="${CONDA_CACHE}/env-${env_id}"
  local donor=""

  if [[ -x "${prefix}/bin/${check_bin}" ]]; then
    echo "    ${yml}: already ready (${prefix})"
    bash "${ROOT}/scripts/fix_easyfuse_pyeasyfuse_env.sh" >/dev/null 2>&1 || true
    return 0
  fi

  # Offline reuse: any existing env-* that already has the required binary.
  for donor in "${CONDA_CACHE}"/env-*; do
    [[ -d "${donor}" ]] || continue
    [[ "${donor}" == "${prefix}" ]] && continue
    if [[ -x "${donor}/bin/${check_bin}" ]]; then
      echo "==> Reusing EasyFuse ${yml}: ${prefix} <= ${donor} (no mamba download)"
      rm -rf "${prefix}"
      ln -sfn "${donor}" "${prefix}"
      if [[ -x "${prefix}/bin/${check_bin}" ]]; then
        bash "${ROOT}/scripts/fix_easyfuse_pyeasyfuse_env.sh" >/dev/null 2>&1 || true
        return 0
      fi
      echo "WARN: symlink reuse failed for ${prefix}; falling through" >&2
      break
    fi
  done

  echo "ERROR: EasyFuse ${yml} needs ${check_bin} at ${prefix}, and no local donor env was found." >&2
  echo "       Refusing mamba env create (intranet cannot reach conda.anaconda.org)." >&2
  echo "       Populate ${CONDA_CACHE} from shared_refs/easyfuse_nextflow_conda or a known-good host cache." >&2
  return 1
}


exec >> "${LOG}" 2>&1
echo ""
echo "==> easyfuse_prebuild_remaining_envs $(date -Is)"

prebuild_conda_env \
  "6f2b394c864eeaa5-8f88fe4572f59d9bb818f7644ca8f1fa" \
  "alignment.yml" \
  "STAR"

prebuild_conda_env \
  "e4b7dd3f8c4b23e0-ae78dd87a2293bf5178c1752ac54f434" \
  "starfusion.yml" \
  "STAR-Fusion"

prebuild_conda_env \
  "cdd3345c1c7c6ebb-12bec81221d755f4cf369bdee252c72d" \
  "fusioncatcher.yml" \
  "fusioncatcher"

prebuild_conda_env \
  "32b8951a86fd0d30-0412febc37d8c32b8c4b8283292cffe1" \
  "requantification.yml" \
  "STAR"

echo "==> easyfuse_prebuild_remaining_envs done $(date -Is)"
