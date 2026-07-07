#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J run_final_TGFM_CD[1-50]
#BSUB -o /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/run_final_TGFM_CD-%I-%J-output.log
#BSUB -e /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/run_final_TGFM_CD-%I-%J-error.log 

# Create Softpack environment
module load HGI/softpack/users/jh59/test-tgfm/2

# Set paths and variables
workdir=/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/output_data/scripts
cd $workdir
TGFM_CODE_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/output_data/TGFM_scripts"
TGFM_INPUT_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/output_data/input_data_for_TGFM_per_chunk"
TiSSUE_SUMMARY_FILE="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/tissue_summary_file.txt"
TISSUE_PRIOR_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/output_data/tissue_prior_output"
OUTPUT_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/output_data/final_TGFM_output"
mkdir -p $OUTPUT_DIR

# Set the jobidentifier as chunk number (loop through each chunk)
CHUNK_NUM=${LSB_JOBINDEX}
echo "Processing chunk $CHUNK_NUM..."

# Run the tissue-specific prior script
python $TGFM_CODE_DIR/run_tgfm_updated.py \
    --trait-name "CD" \
    --tgfm-input-data ${TGFM_INPUT_DIR}/tgfm_chunk${CHUNK_NUM}_input_data_summary.txt \
    --parallel-job-identifier ${CHUNK_NUM} \
    --tissue-summary-file $TiSSUE_SUMMARY_FILE \
    --tissue-specific-prior $TISSUE_PRIOR_DIR/CD \
    --out $OUTPUT_DIR/CD \

echo "Finished processing final TGFM"