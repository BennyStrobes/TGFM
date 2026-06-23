#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J generate_npy_variant_info[1-22]
#BSUB -o /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/generate_npy_variant_info[1-22]-%J-output.log
#BSUB -e /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/generate_npy_variant_info[1-22]-%J-error.log 

# Create Softpack environment
module load HGI/softpack/users/jh59/test-tgfm/1

# Define input/output directories
LD_MATRIX_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file/LD_matrix_file_perwindow"
PLINK_FILE_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file/plink_perwindow"
NPY_OUTPUT_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file/ld_npy_perwindow"              
VARIANT_INFO_OUTPUT_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file/variant_info_file_perwindow"  

# Get the chromosome number from the job index
CHR=${LSB_JOBINDEX}

echo "Processing all windows for $CHR..."

# Loop over all windows belonging to this chromosome
for vcor2_file in ${LD_MATRIX_DIR}/${CHR}:*.unphased.vcor2; do
    if [ -f "$vcor2_file" ]; then
        # Extract window name (remove directory and extension)
        WINDOW_NAME=$(basename "$vcor2_file" .unphased.vcor2)
        
        echo "Processing window: $WINDOW_NAME..."
        python execution_scripts_BH/generate_npy_variant_info.py "$LD_MATRIX_DIR" "$PLINK_FILE_DIR" "$NPY_OUTPUT_DIR" "$VARIANT_INFO_OUTPUT_DIR" "$WINDOW_NAME"
    fi
done

echo "Finished processing all windows for $CHR"
