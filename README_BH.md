## Bradley April 2026
### Re-running TGFM 
Get the paper code
```bash
git clone https://github.com/BennyStrobes/TGFM.git
```

Get relevant submission scripts
```bash
JLDIR=/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis
DIR=execution_scripts_BH
mkdir -p $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/generate_plink_ldmatrix_perwindow.sh $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/generate_plink_ldmatrix_perwindow_for_nan.sh $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/run_generate_npy_variant_info.sh $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/run_generate_final_windows_file.sh $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/generate_npy_variant_info.py $DIR
cp ${JLDIR}/output_data/TGFM_scripts/tgfm.py tgfm_nan_update.py
# REMOVE LINE 23 from this file. This line: mixture_beta_mu = np.clip(mixture_beta_mu, -1.0, 1.0)
```

0. Check the rows of the GWAS files which cannot generate z-score (missing beta or SE)
```bash
module load $uber
Rscript count_missing_z_gwas.r
```

1. Use a single job to regenerate files per chromosome
```bash
bsub -J "Regen_LD[1-22]" -M 100000 -a "memlimit=True" -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" -o logs/generate_plink_ldmatrix_perwindow_for_nan-%I-%J-output.log -e logs/generate_plink_ldmatrix_perwindow_for_nan-%I-%J-error.log -q normal \
"echo '--- Running generate_plink_ldmatrix_perwindow.sh ---' && bash execution_scripts_BH/generate_plink_ldmatrix_perwindow.sh && echo '--- Running generate_plink_ldmatrix_perwindow_for_nan.sh ---' && bash execution_scripts_BH/generate_plink_ldmatrix_perwindow_for_nan.sh && echo '--- Running run_generate_npy_variant_info.sh ---' && bash execution_scripts_BH/run_generate_npy_variant_info.sh"
```

2. Resubmit tgfm script (50 chunks per trait)
```bash
for TRAIT in CD UC IBD; do
    bsub -J "run_final_TGFM_${TRAIT}[1-50]" -M 30000 -a "memlimit=True" -R "select[mem>30000] rusage[mem=30000] span[hosts=1]" -o logs/run_final_TGFM-${TRAIT}-%I-%J-output.log -e logs/run_final_TGFM-${TRAIT}-%I-%J-error.log -q long "bash execution_scripts_BH/run_TGFM-BH-MCP.sh -t ${TRAIT}"
done
```
