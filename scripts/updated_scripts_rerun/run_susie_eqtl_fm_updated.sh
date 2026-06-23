#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J susie_eqtl_fm[1-22]
#BSUB -o /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/susie_eqtl_fm-%I-%J-output.log
#BSUB -e /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/susie_eqtl_fm-%I-%J-error.log 

# Create Softpack environment
module load HGI/softpack/users/jh59/test-tgfm/2

CHROM_NUM=${LSB_JOBINDEX} 

# Set paths and variables
workdir=/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/output_data/scripts
cd $workdir
TGFM_CODE_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/output_data/TGFM_scripts"
PLINK_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/core_analysis_output/IBDverse_multi-tissue_eQTL_project/IBDverse_genotypes/2024_07_11-genotype_plate12345/plink" 
EQTL_SUMSTAT_BASE="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/IBDverse_data/reformatted_eqtl_sumstat" 
GWAS_SUMSTAT_BASE="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/IBDverse_data/DeLange_sumstat" 
OUTPUT_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/output_data/eQTL_fm_chr${CHROM_NUM}" 
mkdir -p "$OUTPUT_DIR"


# list of folders containing eQTL data for all cell type clusters
CELL_TYPE_FOLDERS=(${EQTL_SUMSTAT_BASE}/dMean__*_ct_all)

START_FROM="dMean__T_0_ct_all"


# Loop through each chromosome
echo "Processing chromosome $CHROM_NUM..."

# Flag to indicate when we've found our starting folder
found_start=false

# Loop through each cell type folder
for FOLDER in "${CELL_TYPE_FOLDERS[@]}"; do

    # Get just the folder name without path
    folder_name=$(basename "$FOLDER")
    
    # If we haven't found the starting folder yet, check if this is it
    if [ "$found_start" = false ] && [ "$folder_name" != "$START_FROM" ]; then
        continue  # Skip this iteration
    fi
    
    # We've found our starting folder (or we're past it)
    found_start=true


    if [[ -d "$FOLDER" ]]; then
        # Extract CELL_TYPE
        CELL_TYPE=$(basename "$FOLDER" | tr -s '_' | cut -d'_' -f2)

        # Extract CELL_SUBTYPE if available
        if [[ $(basename "$FOLDER") =~ dMean__${CELL_TYPE}_(.*)_ct_all ]]; then
            CELL_SUBTYPE="${BASH_REMATCH[1]}"
        else
            CELL_SUBTYPE=""
        fi

        # Define eQTL summary statistic file
        if [ -z "$CELL_SUBTYPE" ]; then
            EQTL_FOLDER="${EQTL_SUMSTAT_BASE}/dMean__${CELL_TYPE}_ct_all"
            FILE_PREFIX="${CELL_TYPE}_ct"
        else
            EQTL_FOLDER="${EQTL_SUMSTAT_BASE}/dMean__${CELL_TYPE}_${CELL_SUBTYPE}_ct_all"
            FILE_PREFIX="${CELL_TYPE}_${CELL_SUBTYPE}_ct"
        fi
    
        EQTL_SUMSTAT="${EQTL_FOLDER}/reformatted_cis_nominal1_eqtl.${CHROM_NUM}.tsv"
        #GWAS_SUMSTAT="${GWAS_SUMSTAT_BASE}/InflammatoryBowelDisease_DeLange_NatGen2017_formatted.txt.gz"
        GWAS_SUMSTAT="${GWAS_SUMSTAT_BASE}/IBD_geno_intersected_GWAS_sumstat.txt.gz"
        QVAL_FILE="${EQTL_FOLDER}/Cis_eqtls_qval.tsv"

        EQTL_OUTPUT_STEM="${OUTPUT_DIR}/${FILE_PREFIX}"

        # Define PLINK genotype file stem
        PLINK_GENO_FILE_STEM="${PLINK_DIR}/imputed_chr${CHROM_NUM}"

        # Run the fine-mapping script
        python "${TGFM_CODE_DIR}/susie_eqtl_fine_mapping_for_tgfm.py" \
            --eqtl-data-type SumStat \
            --chrom "$CHROM_NUM" \
            --genotype-stem "$PLINK_GENO_FILE_STEM" \
            --eqtl-sumstat "$EQTL_SUMSTAT" \
            --gwas-sumstat "$GWAS_SUMSTAT" \
            --filter-strand-ambiguous \
            --out "$EQTL_OUTPUT_STEM" \
            --qval-file "$QVAL_FILE" \
                
        echo "Completed: Chromosome $CHROM_NUM, CellType: $CELL_TYPE, CellSubtype: $CELL_SUBTYPE"
        
    fi
done


echo "All tasks completed!"