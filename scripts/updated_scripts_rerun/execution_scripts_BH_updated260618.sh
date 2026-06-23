#Re-running TGFM
#Get the paper code
git clone https://github.com/BennyStrobes/TGFM.git

#Get relevant submission scripts
JLDIR=/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis
DIR=execution_scripts_BH
mkdir -p $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/generate_plink_ldmatrix_perwindow.sh $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/generate_plink_ldmatrix_perwindow_for_nan.sh $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/run_generate_npy_variant_info.sh $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/run_generate_final_windows_file.sh $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/generate_npy_variant_info.py $DIR
cp ${JLDIR}/output_data/TGFM_scripts/tgfm.py tgfm_nan_update_v2.py 
# REMOVE LINE 23 from this file. This line: mixture_beta_mu = np.clip(mixture_beta_mu, -1.0, 1.0)

cp ${JLDIR}/prepare_input_data/scripts/filter_GWAS_genotype/filter_GWAS_genotype_intersected.sh $DIR
cp ${JLDIR}/prepare_input_data/scripts/filter_GWAS_genotype/correct_flip.sh $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/generate_windows.py $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/run_generate_windows.sh $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/make_nan_indices.py $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/make_nan_indices.sh $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/generate_final_windows_file.py $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/run_generate_final_windows_file.sh $DIR

cp ${JLDIR}/output_data/TGFM_scripts/susie_eqtl_fine_mapping_for_tgfm.py $DIR
cp ${JLDIR}/output_data/make_tissue_summary_file.py $DIR
cp ${JLDIR}/output_data/TGFM_scripts/run_tgfm_without_sampling.py $DIR
cp ${JLDIR}/output_data/TGFM_scripts/run_tgfm_tissue_specific_prior.py $DIR
cp ${JLDIR}/output_data/run_tissue_prior.sh $DIR

# intersect the GWAS genotype and GWAS sumstat (intersected of three traits) first
bsub < execution_scripts_BH/intersect_GWAS_geno_sumsta.sh

# re-run the GWAS genotype data filtering according to the intersected gwas sumstat
rm /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/output_data/scripts/check_flip_output/check_flip_output/*.txt
bash execution_scripts_BH/correct_flip.sh


#Use a single job to regenerate files per chromosome (starting from generating windows, to generate ld, generate nan list, re-run the ld generation for nan, finally merge the npy and variant info to the window files)
bsub -J "Regen_LD[1-22]" \
     -M 100000 \
     -a "memlimit=True" \
     -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" \
     -o logs/generate_plink_ldmatrix_perwindow_for_nan-%I-%J-output.log \
     -e logs/generate_plink_ldmatrix_perwindow_for_nan-%I-%J-error.log \
     -q normal \
     "
echo '--- Running run_generate_windows.sh ---' && \
bash exexecution_scripts_BH/run_generate_windows.sh && \

echo '--- Running generate_plink_ldmatrix_perwindow.sh ---' && \
bash execution_scripts_BH/generate_plink_ldmatrix_perwindow.sh && \

echo '--- Running make_nan_indices.sh ---' && \
bash execution_scripts_BH/make_nan_indices.sh && \

echo '--- Running generate_plink_ldmatrix_perwindow_for_nan.sh ---' && \
bash execution_scripts_BH/generate_plink_ldmatrix_perwindow_for_nan.sh && \

echo '--- Running run_generate_npy_variant_info.sh ---' && \
bash execution_scripts_BH/run_generate_npy_variant_info.sh && \

echo '--- Running run_generate_final_windows_file.sh ---' && \
bash execution_scripts_BH/run_generate_final_windows_file.sh.sh
"

# run TGFM 

# Re-run the eQTL fine-mapping (filtering eQTL variants by intersected GWAS sumstat)
bsub < execution_scripts_BH/run_susie_eqtl_fm_updated.sh

# re-run the data preparation (integration with the updated eQTL fine-mapping result) + randomly sample the updated eQTL fine-mapping for each gt pair
python execution_scripts_BH/make_tissue_summary_file.py
python execution_scripts_BH/make_gwas_summary_file.py
bsub < execution_scripts_BH/run_prepare_input_updated.sh

# re-run the simple version of TGFM (without random sampling)
bsub < execution_scripts_BH/generate_window_per_chunk.sh # make 50 chunks
bsub < execution_scripts_BH/run_TGFM_without_sampling.sh # repeat it for all three traits (IBD, UC, CD), but recommend to do a test on one chunk first

# re-run the prior distribution based on the result from TGFM without sampling
bsub < execution_scripts_BH/run_tissue_prior.sh # repeat it for all three traits

# re-run the final TGFM (with 100 randomly sampled expression model from eQTL fine-mapping)
bsub < execution_scripts_BH/run_TGFM.sh # repeat it for all three traits

#Resubmit tgfm script (50 chunks per trait)
for TRAIT in CD UC IBD; do
    bsub -J "run_final_TGFM_${TRAIT}[1-50]" -M 30000 -a "memlimit=True" -R "select[mem>30000] rusage[mem=30000] span[hosts=1]" -o logs/run_final_TGFM-${TRAIT}-%I-%J-output.log -e logs/run_final_TGFM-${TRAIT}-%I-%J-error.log -q long "bash execution_scripts_BH/run_TGFM-BH-MCP.sh -t ${TRAIT}"
done

