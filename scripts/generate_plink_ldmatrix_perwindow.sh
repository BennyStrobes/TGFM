#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J generate_plink_ldmatrix_perwindow[10]
#BSUB -o /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/generate_plink_ldmatrix_perwindow-%I-%J-output.log
#BSUB -e /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/generate_plink_ldmatrix_perwindow-%I-%J-error.log

# Set working directory and directory to latest version of PLINK2
workdir=/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/scripts/generate_window_file
cd $workdir
plink2_bin=/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/scripts/generate_window_file/plink2
export PATH="$workdir:$PATH"
chmod +x plink2

# Run PLINK2
$plink2_bin --version

# Set paths and variables
plink_prefix="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/human_core_exome_genotype/humancoreexome_filtered_plink_corrected_flip"  
window_file="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file/windows/windows_chr${LSB_JOBINDEX}.tsv"
plink_output_dir="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file/plink_perwindow"
ld_output_dir="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file/LD_matrix_file_perwindow"

# List of chromosomes to loop through (defined by job index)
CHR=${LSB_JOBINDEX}

# Start from a specific window
#START_WINDOW="2:235000001-238000001"

# Flag to control when to start processing
#start_processing=false

while IFS=$'\t' read -r window_name chr start end; do
    # Skip header line
    [[ "$window_name" == "window_name" ]] && continue

    # Check if current window matches the START_WINDOW
    #if [[ "$window_name" == "$START_WINDOW" ]]; then
    #    start_processing=true
    #fi

    # Skip lines until the START_WINDOW is found
    #if [[ "$start_processing" == false ]]; then
    #    continue
    #fi
    
    # Extract PLINK data for the window
    plink2 --bfile "$plink_prefix" \
          --chr "$chr" \
          --from-bp "$start" \
          --to-bp "$end" \
          --make-bed \
          --out "$plink_output_dir/$window_name"

    # Compute LD matrix
    plink2 --bfile "$plink_output_dir/$window_name" \
           --r2-unphased square 'ref-based' \
           --out "$ld_output_dir/$window_name"

done < "$window_file"
