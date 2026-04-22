## Bradley April 2026
### Re-running TGFM 
Get the paper code
```bash
git clone https://github.com/BennyStrobes/TGFM.git
```

Get relevant submission scripts
```
JLDIR=/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis
DIR=execution_scripts_BH
mkdir -p $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/generate_plink_ldmatrix_perwindow.sh $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/generate_plink_ldmatrix_perwindow_for_nan.sh $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/run_generate_npy_variant_info.sh $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/run_generate_final_windows_file.sh $DIR
cp ${JLDIR}/prepare_input_data/scripts/generate_window_file/generate_npy_variant_info.py $DIR
```

1. Use a single job to regenerate files per chromosome
```bash
bsub -J "Regen_LD[19]" -M 100000 -a "memlimit=True" -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" -o logs/generate_plink_ldmatrix_perwindow_for_nan-%I-%J-output.log -e logs/generate_plink_ldmatrix_perwindow_for_nan-%I-%J-error.log -q normal \
"echo '--- Running generate_plink_ldmatrix_perwindow.sh ---' && bash execution_scripts_BH/generate_plink_ldmatrix_perwindow.sh && echo '--- Running generate_plink_ldmatrix_perwindow_for_nan.sh ---' && bash execution_scripts_BH/generate_plink_ldmatrix_perwindow_for_nan.sh && echo '--- Running run_generate_npy_variant_info.sh ---' && bash execution_scripts_BH/run_generate_npy_variant_info.sh"
```

1. Regenerate plink per window
```
bsub -J "generate_plink_ldmatrix_perwindow_for_nan[22]" -M 10000 -a "memlimit=True" -R "select[mem>10000] rusage[mem=10000] span[hosts=1]" -o logs/generate_plink_ldmatrix_perwindow_for_nan-%I-%J-output.log -e logs/generate_plink_ldmatrix_perwindow_for_nan-%I-%J-error.log -q normal "bash execution_scripts_BH/generate_plink_ldmatrix_perwindow_for_nan.sh"
```

2. Regenerate plink per window for nan
```
bsub -J "generate_plink_ldmatrix_perwindow_for_nan[22]" -M 10000 -a "memlimit=True" -R "select[mem>10000] rusage[mem=10000] span[hosts=1]" -o logs/generate_plink_ldmatrix_perwindow_for_nan-%I-%J-output.log -e logs/generate_plink_ldmatrix_perwindow_for_nan-%I-%J-error.log -q normal "bash execution_scripts_BH/generate_plink_ldmatrix_perwindow_for_nan.sh"
```

3. Convert vcor to npy files
```
bsub -J "run_generate_npy_variant_info[1-22]" -M 100000 -a "memlimit=True" -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" -o logs/run_generate_npy_variant_info-%I-%J-output.log -e logs/run_generate_npy_variant_info-%I-%J-error.log -q normal "bash execution_scripts_BH/run_generate_npy_variant_info.sh"
```

2. Resubmit tgfm script (1-50 chunks, 1 for each trait)
```
TRAIT=CD
bsub -J "run_final_TGFM[1]" -M 100000 -a "memlimit=True" -R "select[mem>100000] rusage[mem=100000] span[hosts=1]" -o logs/run_final_TGFM-${TRAIT}-%I-%J-output.log -e logs/run_final_TGFM-${TRAIT}-%I-%J-error.log -q normal "bash execution_scripts_BH/run_TGFM-BH-MCP.sh -t ${TRAIT}"

```
