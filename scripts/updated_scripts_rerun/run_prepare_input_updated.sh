#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J run_prepare_input[1-22]
#BSUB -o /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/run_prepare_input-%I-%J-output.log
#BSUB -e /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/run_prepare_input-%I-%J-error.log 

# Create Softpack environment
module load HGI/softpack/users/jh59/test-tgfm/2

# Set paths and variables
workdir=/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/output_data/scripts
cd $workdir
TGFM_CODE_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/output_data/TGFM_scripts"
WINDOW_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file"
TiSSUE_SUMMARY_FILE="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/tissue_summary_file.txt"
GWAS_SUMMARY_FILE="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/gwas_summary_file.txt"
OUTPUT_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/output_data/input_data_for_TGFM"

CHROM=${LSB_JOBINDEX}
echo "Processing chromosome $CHROM..."

# Run the prepare_input.py script
python ${TGFM_CODE_DIR}/prepare_input_data_for_tgfm_updated.py \
    --chrom ${CHROM} \
    --window-file ${WINDOW_DIR}/window_file_chr${CHROM}.tsv \
    --tissue-summary-file ${TiSSUE_SUMMARY_FILE} \
    --gwas-summary-file ${GWAS_SUMMARY_FILE} \
    --standardize-gwas-summary-statistics \
    --cis-window-size 1000000 \
    --out ${OUTPUT_DIR}/tgfm

echo "Finished processing chromosome $CHROM"