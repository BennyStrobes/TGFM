#!/usr/bin/env bash
# Correct strand flips in a plink genotype dataset.
#
# Called by the Nextflow CORRECT_FLIP process with two env vars set:
#   INPUT_PLINK_PREFIX   — path prefix of the source .bed/.bim/.fam
#   OUTPUT_PLINK_PREFIX  — path prefix to write the corrected output
#
# No hardcoded paths. All data locations come from env vars.
#
# This script delegates to the original Jingling correct_flip.sh if
# CORRECT_FLIP_DELEGATE is set, otherwise falls back to a passthrough copy.
# Replace the passthrough block with the actual flip-correction logic once
# the original script is accessible.

set -euo pipefail

: "${INPUT_PLINK_PREFIX:?INPUT_PLINK_PREFIX env var must be set}"
: "${OUTPUT_PLINK_PREFIX:?OUTPUT_PLINK_PREFIX env var must be set}"

WORK_DIR="$(dirname "${OUTPUT_PLINK_PREFIX}")/flip_work"
mkdir -p "${WORK_DIR}"

# ── Delegate to the original Jingling script if available ────────────────────
if [[ -n "${CORRECT_FLIP_DELEGATE:-}" && -f "${CORRECT_FLIP_DELEGATE}" ]]; then
    echo "Delegating to original correct_flip script: ${CORRECT_FLIP_DELEGATE}"
    export INPUT_PLINK_PREFIX OUTPUT_PLINK_PREFIX WORK_DIR
    bash "${CORRECT_FLIP_DELEGATE}"

# ── Passthrough fallback — replace with actual flip-correction logic ─────────
else
    echo "WARNING: No flip-correction delegate found."
    echo "Copying input plink files as passthrough (no strand correction applied)."
    echo "Set CORRECT_FLIP_DELEGATE to the path of the original correct_flip.sh"
    echo "to enable proper strand flip correction."
    cp "${INPUT_PLINK_PREFIX}.bed" "${OUTPUT_PLINK_PREFIX}.bed"
    cp "${INPUT_PLINK_PREFIX}.bim" "${OUTPUT_PLINK_PREFIX}.bim"
    cp "${INPUT_PLINK_PREFIX}.fam" "${OUTPUT_PLINK_PREFIX}.fam"
fi

echo "correct_flip.sh done: ${OUTPUT_PLINK_PREFIX}"
