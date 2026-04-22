#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J generate_final_windows_file[1-22]
#BSUB -o /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/generate_final_windows_file[1-22]-%J-output.log
#BSUB -e /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/generate_final_windows_file[1-22]-%J-error.log 

# Create Softpack environment
module load HGI/softpack/users/jh59/test-tgfm/1

# Set working directory
workdir=/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/scripts/generate_window_file
cd $workdir

# Define input/output directories
window_input_dir=
WINDOW_INPUT_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file/windows"
NPY_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file/ld_npy_perwindow"                  
VARIANT_INFO_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file/variant_info_file_perwindow"   
OUTPUT_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file"        

# List of chromosomes to loop through
CHROMOSOME=${LSB_JOBINDEX} 
   
echo "Processing chromosome $CHROMOSOM..."
python generate_final_windows_file.py "$CHROMOSOME" "$WINDOW_INPUT_DIR" "$NPY_DIR" "$VARIANT_INFO_DIR" "$OUTPUT_DIR" 

echo "Chromosome $CHROMOSOME processed! Final window files stored in $OUTPUT_DIR."