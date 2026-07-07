#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J generate_window_per_chunk
#BSUB -o /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/generate_window_per_chunk-%J-output.log
#BSUB -e /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/generate_window_per_chunk-%J-error.log 

# Create Softpack environment
module load HGI/softpack/users/jh59/test-tgfm/2


# Set paths and variables
workdir=/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/output_data/scripts/generate_window_per_chunk
cd $workdir
CODE_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/output_data/scripts/generate_window_per_chunk"
INPUT_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/output_data/input_data_for_TGFM"
OUTPUT_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/output_data/input_data_for_TGFM_per_chunk"
mkdir -p ${OUTPUT_DIR}

# Run the generate_window_per_chunk.py script
python ${CODE_DIR}/generate_window_per_chunk.py

echo "Finished processing window file per chunk"