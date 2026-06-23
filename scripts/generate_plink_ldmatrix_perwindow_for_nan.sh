#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J generate_plink_ldmatrix_perwindow_for_nan[1-22]
#BSUB -o /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/generate_plink_ldmatrix_for_nan[1-22]-%J-output.log
#BSUB -e /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/generate_plink_ldmatrix_for_nan[1-22]-%J-error.log

# Set working directory and directory to latest version of PLINK2
workdir=/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/scripts/generate_window_file
cd $workdir
plink2_bin=/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/scripts/generate_window_file/plink2
export PATH="$workdir:$PATH"
#chmod +x plink2

# Run PLINK2
$plink2_bin --version

# Set paths and variables
exclude_file="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file/LD_matrix_file_perwindow/nan_indices/total_nan_indices.txt"
plink_prefix="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/human_core_exome_genotype/humancoreexome_filtered_plink_corrected_flip"  
window_file="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file/windows/windows_chr${LSB_JOBINDEX}.tsv"
plink_output_dir="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file/plink_perwindow"
ld_output_dir="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/prepare_input_data/window_file/LD_matrix_file_perwindow"

# List of chromosomes to loop through (defined by job index)
CHR=${LSB_JOBINDEX}

# extract the windows contains nan values in LD matrix
declare -A WINDOW_WITH_NAN
while IFS=$'\t' read -r window variant; do
    WINDOW_WITH_NAN["$window"]="${WINDOW_WITH_NAN["$window"]} $variant"
done < "$exclude_file"


while IFS=$'\t' read -r window_name chr start end; do
    # Skip header line
    [[ "$window_name" == "window_name" ]] && continue

     # Check if this window contain variants with nan in LD matrix
    if [[ -z "${WINDOW_WITH_NAN[$window_name]}" ]]; then
        continue  # if it doesn't contain nan, Skip this window
    fi

    echo "Processing window: $window_name"
    # Create a temporary exclude file for this window
    exclude_temp_file=$(mktemp)
    for v in ${WINDOW_WITH_NAN["$window_name"]}; do
        echo "$v"
    done > "$exclude_temp_file"

    # Extract PLINK data for the window
    plink2 --bfile "$plink_prefix" \
          --chr "$chr" \
          --from-bp "$start" \
          --to-bp "$end" \
          --exclude "$exclude_temp_file" \
          --make-bed \
          --out "$plink_output_dir/$window_name"

    # Compute LD matrix
    plink2 --bfile "$plink_output_dir/$window_name" \
           --r2-unphased square 'ref-based' \
           --out "$ld_output_dir/$window_name"

    # Cleanup
    rm "$exclude_temp_file"

done < "$window_file"
