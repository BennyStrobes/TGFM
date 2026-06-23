# TGFM Pipeline Overview

Tissue-specific Gene Fine-Mapping (TGFM) for IBD traits (CD, UC, IBD) using eQTL fine-mapping + GWAS summary statistics. All outputs land in `results/`.

---

## Stage 1 — Genotype & GWAS Preparation

### INTERSECT_GWAS_GENO
Finds the set of SNPs common to all three GWAS traits (CD, UC, IBD) *and* present in the human core exome genotype array. Filters both the genotype plink files and each GWAS summary statistic file down to this shared SNP set.

- **Input**: raw humancore exome plink, three GWAS sumstat `.gz` files, pre-computed common-SNP list
- **Output**: `results/genotype/humancoreexome_filtered_plink.{bed,bim,fam}`, `results/sumstats/{IBD,CD,UC}_geno_intersected_GWAS_sumstat.txt.gz`

### CORRECT_FLIP
Corrects strand flips in the filtered genotype (A/T and C/G ambiguous SNPs that could be on either strand). Calls Jingling's `correct_flip.sh`.

- **Input**: `results/genotype/humancoreexome_filtered_plink.*`
- **Output**: `results/genotype/humancoreexome_filtered_plink_corrected_flip.{bed,bim,fam}`

---

## Stage 2 — LD Matrix Generation (per chromosome, run in parallel for chr 1–22)

### GENERATE_WINDOWS
Divides each chromosome into overlapping 3 Mb windows (with 1 Mb step) anchored to the genotype variants. Produces the window definition file used by all downstream LD steps.

- **Input**: corrected-flip plink, chromosome number
- **Output**: `results/window_files/windows_chr{N}.tsv`

### GENERATE_LD_MATRIX
For each window on a chromosome, extracts a per-window plink bed and computes the LD (r²) matrix using plink2. Windows that fall in centromeric gaps (no variants) are silently skipped rather than failing the whole job.

- **Input**: corrected-flip plink, window TSV
- **Output**: `results/ld_prep/plink_perwindow/{window}.{bed,bim,fam}`, `results/ld_prep/LD_matrix_perwindow/{window}.unphased.vcor2`

### MAKE_NAN_INDICES
Scans the raw LD matrices for NaN entries (caused by monomorphic or near-monomorphic variants) and records which variant in which window caused each NaN.

- **Input**: `results/ld_prep/LD_matrix_perwindow/`, chromosome number
- **Output**: `results/ld_prep/nan_indices/nan_indices_chr{N}.txt`

### COLLECT_NAN_INDICES
Concatenates all per-chromosome NaN index files into a single lookup file.

- **Input**: all `nan_indices_chr*.txt` files
- **Output**: `results/ld_prep/nan_indices/total_nan_indices.txt`

### FIX_NAN_LD_MATRIX
Re-computes LD for any window that had NaN entries, this time excluding the offending variants. Produces clean LD matrices.

- **Input**: corrected-flip plink, window TSV, `total_nan_indices.txt`
- **Output**: updated files in `results/ld_prep/plink_perwindow/` and `results/ld_prep/LD_matrix_perwindow/`

### GENERATE_NPY_VARIANT_INFO
Converts the text-format `.vcor2` LD matrices to compressed NumPy `.npy` arrays and extracts per-variant metadata (chromosome, position, alleles) for each window.

- **Input**: `results/ld_prep/LD_matrix_perwindow/`, `results/ld_prep/plink_perwindow/`
- **Output**: `results/ld_prep/ld_npy_perwindow/{window}_ld.npy`, `results/ld_prep/variant_info_perwindow/{window}_variant_info.txt`

### GENERATE_FINAL_WINDOWS
Merges the window definitions with the paths to the `.npy` LD files and variant info files into a single final window file used by the TGFM input preparation step.

- **Input**: window TSV, npy dir, variant info dir
- **Output**: `results/window_files/window_file_chr{N}.tsv`

---

## Stage 3 — eQTL Fine-Mapping (per chromosome, run in parallel)

### SUSIE_EQTL_FM
Runs SuSiE eQTL fine-mapping for every cell type separately, on each chromosome. Uses the intersected GWAS sumstat to filter eQTL variants to the same SNP set as the GWAS. Produces posterior inclusion probabilities (PIPs) and effect-size distributions for each gene–tissue pair.

- **Input**: per-cell-type eQTL sumstats (`dMean__*_ct_all/`), eQTL imputed genotype, IBD GWAS sumstat
- **Script**: `scripts/susie_eqtl_fine_mapping_for_tgfm.py`
- **Output**: `results/eqtl_fine_mapping/{cell_type}_ct_chr{N}_gene_summary.txt` + associated `.npy` files (one set per cell type × chromosome)

---

## Stage 4 — Summary File Construction

### MAKE_TISSUE_SUMMARY
Scans the eQTL fine-mapping outputs to discover which cell types completed, and writes a lookup table mapping each tissue/cell type to its sample size and the file path stem used by TGFM.

- **Input**: `results/eqtl_fine_mapping/`
- **Output**: `results/tissue_summary_file.txt` (`tissue_name`, `sample_size`, `output_stem` columns)

### MAKE_GWAS_SUMMARY
Reads the intersected GWAS sumstat files to extract sample sizes and write a lookup table mapping each trait to its sumstat file and sample size.

- **Input**: `results/sumstats/`
- **Script**: `scripts/make_gwas_summary_file.py`
- **Output**: `results/gwas_summary_file.txt` (`trait_name`, `sample_size`, `summary_statistics` columns)

---

## Stage 5 — TGFM Input Preparation

### PREPARE_INPUT
For each chromosome, integrates the final window file, eQTL fine-mapping posteriors, and GWAS summary statistics into a single structured input object that TGFM can consume. Performs allele alignment and optionally standardises GWAS effect sizes.

- **Input**: `window_file_chr{N}.tsv`, `tissue_summary_file.txt`, `gwas_summary_file.txt`, LD npy files
- **Script**: `scripts/prepare_input_data_for_tgfm_updated.py`
- **Output**: `results/tgfm_input/tgfm_chr{N}_input_data_summary.txt` + associated data files

### GENERATE_WINDOW_PER_CHUNK
Sorts all windows genome-wide and splits them into 50 equal-sized chunks for parallel downstream processing.

- **Input**: `results/tgfm_input/tgfm_chr*_input_data_summary.txt`
- **Script**: `scripts/generate_window_per_chunk.py`
- **Output**: `results/tgfm_input_chunks/tgfm_chunk{1..50}_input_data_summary.txt`

---

## Stage 6 — TGFM Model Fitting (per trait × chunk)

### TGFM_WITHOUT_SAMPLING
Runs a fast version of TGFM (no random sampling of eQTL effect sizes) on each of the 50 chunks × 3 traits (150 jobs). Used to learn the tissue-specific prior distribution over gene–tissue effect sizes.

- **Input**: chunk summary file, GWAS trait
- **Script**: `scripts/run_tgfm_without_sampling.py`
- **Output**: `results/tgfm_without_sampling/{trait}/` — per-window TGFM results without sampling

### TISSUE_PRIOR
Aggregates the no-sampling TGFM results across all 50 chunks for a given trait to estimate the prior probability that each tissue is causal. Run once per trait after all 50 chunks complete.

- **Input**: `results/tgfm_without_sampling/{trait}/`, chunk ID list, tissue summary file
- **Script**: `scripts/run_tgfm_tissue_specific_prior.py`
- **Output**: `results/tissue_prior/{trait}/` — fitted prior distribution over tissues

### FINAL_TGFM
Runs the full TGFM model with 100 random samples of eQTL expression effects, using the tissue-specific prior from the previous step. Produces final gene–tissue posterior inclusion probabilities (PIPs) for each window. 50 chunks × 3 traits = 150 jobs.

- **Input**: chunk summary file, tissue summary file, `results/tissue_prior/{trait}/`
- **Script**: `scripts/run_tgfm_updated.py` (imports `tgfm_nan_update_v2.py`)
- **Output**: `results/tgfm_final/{trait}/` — final PIPs and credible sets per window

---

## Output Directory Summary

```
results/
├── genotype/            # Filtered & strand-corrected genotype plink files
├── sumstats/            # GWAS sumstats filtered to the common SNP set
├── window_files/        # Window definitions (raw + final with LD paths)
├── ld_prep/
│   ├── plink_perwindow/         # Per-window plink beds
│   ├── LD_matrix_perwindow/     # Raw LD matrices (.vcor2)
│   ├── ld_npy_perwindow/        # Compressed LD arrays (.npy)
│   ├── variant_info_perwindow/  # Per-window variant metadata
│   └── nan_indices/             # NaN variant lookup tables
├── eqtl_fine_mapping/   # SuSiE posteriors per cell type × chromosome
├── tgfm_input/          # Per-chromosome TGFM input objects
├── tgfm_input_chunks/   # Above split into 50 genome-wide chunks
├── tgfm_without_sampling/{CD,UC,IBD}/  # Fast TGFM for prior estimation
├── tissue_prior/{CD,UC,IBD}/           # Fitted tissue-specific priors
├── tgfm_final/{CD,UC,IBD}/             # Final gene–tissue PIPs (main result)
├── tissue_summary_file.txt
├── gwas_summary_file.txt
└── pipeline_info/       # Nextflow execution reports and traces
```
