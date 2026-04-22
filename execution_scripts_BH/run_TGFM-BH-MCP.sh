#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J run_final_TGFM[16]
#BSUB -o /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/run_final_TGFM-%I-%J-output.log
#BSUB -e /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/run_final_TGFM-%I-%J-error.log

usage() {
    echo "Usage: $0 -t TRAIT"
    echo "  -t  Trait name (e.g. CD)  (required)"
    echo "  -h  Show this help message"
    exit 1
}

while getopts ":t:h" opt; do
    case $opt in
        t) TRAIT="$OPTARG" ;;
        h) usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
        \?) echo "Unknown option: -$OPTARG" >&2; usage ;;
    esac
done

if [ -z "${TRAIT}" ]; then
    echo "Error: -t TRAIT is required." >&2
    usage
fi

# Create Softpack environment
module load HGI/softpack/users/jh59/test-tgfm/2

# Set paths and variables
workdir=/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/output_data/scripts # Same as jingling
cd $workdir
TGFM_CODE_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/TGFM" # Original https://github.com/BennyStrobes/TGFM
TGFM_INPUT_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/output_data/input_data_for_TGFM_per_chunk_mcp" # Major cell pop
TiSSUE_SUMMARY_FILE="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/mcp_summary_file.txt" # major cell pop
TISSUE_PRIOR_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/output_data/mcp_prior_output" # major cell pop
OUTPUT_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/TGFM/final_TGFM_output_mcp-BH" # new output dir
mkdir -p $OUTPUT_DIR

# Set the jobidentifier as chunk number (loop through each chunk)
CHUNK_NUM=${LSB_JOBINDEX}
echo "Processing chunk $CHUNK_NUM..."

# Run the tissue-specific prior script
python $TGFM_CODE_DIR/run_tgfm.py \
    --trait-name "${TRAIT}" \
    --tgfm-input-data ${TGFM_INPUT_DIR}/tgfm_chunk${CHUNK_NUM}_input_data_summary.txt \
    --parallel-job-identifier ${CHUNK_NUM} \
    --tissue-summary-file $TiSSUE_SUMMARY_FILE \
    --tissue-specific-prior $TISSUE_PRIOR_DIR/${TRAIT} \
    --out $OUTPUT_DIR/${TRAIT} \

echo "Finished processing final TGFM"

# Submission
# TRAIT=CD
# bsub -J "run_final_TGFM[1-2]" -M 100000 -a "memlimit=True" -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" -o /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/TGFM/logs/run_final_TGFM-${TRAIT}-%I-%J-output.log -e /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/TGFM/logs/run_final_TGFM-${TRAIT}-%I-%J-error.log -q normal "bash run_TGFM-BH-MCP.sh -t ${TRAIT}"

