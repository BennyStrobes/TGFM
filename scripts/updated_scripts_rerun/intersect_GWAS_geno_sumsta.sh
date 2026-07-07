#!/bin/bash
#BSUB -q normal
#BSUB -G team152
#BSUB -n 1
#BSUB -M 100000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -J intersect_GWAS_genotype_sumstat
#BSUB -o /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/intersect_GWAS_genotype_sumstat%J-output.log
#BSUB -e /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/farm_log/intersect_GWAS_genotype_sumstat-%J-error.log 

# Create Softpack environment
module load HGI/softpack/users/jh59/test-tgfm/2

# Set paths and variables
workdir=/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis
cd $workdir

# INPUT
GENO_PREFIX="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/human_exome_core_genotype/humancoreexome_allchr_subset_included_in_ibd_analysis"
SUMSTAT_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/IBDverse_data/DeLange_sumstat"
THREE_TRAIT_SNPS="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis//prepare_input_data/scripts/intersect_GWAS_sumstat/common_snps.txt"
SNP_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis//prepare_input_data/scripts/intersect_GWAS_sumstat"
OUTPUT_DIR="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/human_exome_core_genotype" 
mkdir -p ${OUTDIR}

FILTERED_PLINK_PREFIX="${OUTPUT_DIR}/humancoreexome_filtered_plink"

# Create final SNP list:
# (IBD ∩ CD ∩ UC) ∩ Genotype

cut -f2 ${GENO_PREFIX}.bim | sort > ${SNP_DIR}/geno_snps.txt
sort ${THREE_TRAIT_SNPS} > ${SNP_DIR}/trait_snps.txt

comm -12 \
    ${SNP_DIR}/trait_snps.txt \
    ${SNP_DIR}/geno_snps.txt \
    > ${SNP_DIR}/final_common_snps.txt

echo "Final SNP count:"
wc -l ${SNP_DIR}/final_common_snps.txt


# Filter genotype
plink \
    --bfile ${GENO_PREFIX} \
    --extract ${SNP_DIR}/final_common_snps.txt \
    --make-bed \
    --out ${FILTERED_PLINK_PREFIX}


# Filter all three sumstats
for trait in IBD CD UC
do
    zcat ${SUMSTAT_DIR}/${trait}_intersected_GWAS_sumstat.txt.gz | \
    awk 'BEGIN{
             while((getline < snpfile)>0) keep[$1]=1
         }
         NR==1 || ($2 in keep)' \
         snpfile=${SNP_DIR}/final_common_snps.txt \
    > ${SUMSTAT_DIR}/${trait}_geno_intersected_GWAS_sumstat.txt
done

gzip -f ${SUMSTAT_DIR}/${trait}_geno_intersected_GWAS_sumstat.txt

echo "Done."