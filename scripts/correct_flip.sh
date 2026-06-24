#!/usr/bin/env bash
# Correct strand flips in a PLINK genotype dataset using plink2 --update-alleles.
#
# Required env vars:
#   INPUT_PLINK_PREFIX   — source plink prefix (.bed/.bim/.fam)
#   OUTPUT_PLINK_PREFIX  — destination prefix for corrected output
#   PLINK2_BIN           — path to the plink2 executable
#
# Optional env var:
#   FLIP_ALLELES_FILE    — path to the allele-update file (e.g. combined_flips2.txt).
#                          If not set or the file does not exist, the input plink
#                          files are copied as-is (no strand correction applied).

set -euo pipefail

: "${INPUT_PLINK_PREFIX:?INPUT_PLINK_PREFIX env var must be set}"
: "${OUTPUT_PLINK_PREFIX:?OUTPUT_PLINK_PREFIX env var must be set}"
: "${PLINK2_BIN:?PLINK2_BIN env var must be set}"

mkdir -p "$(dirname "${OUTPUT_PLINK_PREFIX}")"

if [[ -n "${FLIP_ALLELES_FILE:-}" && -f "${FLIP_ALLELES_FILE}" ]]; then
    echo "Correcting strand flips using: ${FLIP_ALLELES_FILE}"
    "${PLINK2_BIN}" \
        --bfile  "${INPUT_PLINK_PREFIX}" \
        --update-alleles "${FLIP_ALLELES_FILE}" \
        --make-bed \
        --out    "${OUTPUT_PLINK_PREFIX}"
    echo "Corrected flip output: ${OUTPUT_PLINK_PREFIX}.bed/bim/fam"
else
    echo "WARNING: FLIP_ALLELES_FILE not set or not found — copying input files without strand correction."
    echo "  Set FLIP_ALLELES_FILE to the combined flips file (e.g. combined_flips2.txt) to enable correction."
    cp "${INPUT_PLINK_PREFIX}.bed" "${OUTPUT_PLINK_PREFIX}.bed"
    cp "${INPUT_PLINK_PREFIX}.bim" "${OUTPUT_PLINK_PREFIX}.bim"
    cp "${INPUT_PLINK_PREFIX}.fam" "${OUTPUT_PLINK_PREFIX}.fam"
fi

echo "correct_flip.sh done: ${OUTPUT_PLINK_PREFIX}"
