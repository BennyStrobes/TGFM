# Tissue-Gene Fine-Mapping (TGFM) — IBDverse

This repository applies TGFM to IBD traits (CD, UC, IBD) using multi-tissue eQTL data from the IBDverse project. It is based on the method published in:

[Strober et al. Fine-mapping causal tissues and genes at disease-associated loci. *Nature Genetics*, 2025.](https://www.nature.com/articles/s41588-024-01994-2)

---

## Changes relative to the original Jingling analysis

The original analysis (`jingling_analysis/`) was run as a series of manual `bsub` array jobs with hardcoded paths. This fork makes the following changes:

### 1. Nextflow pipeline
All steps are re-implemented as a single Nextflow (DSL2) pipeline (`workflow/main.nf`), replacing the collection of individual bsub submission scripts. The pipeline:
- Runs on the Sanger LSF cluster via the `sanger` profile
- Uses Singularity containers throughout (no `module load`)
- Supports `-resume` to restart from any failed step without re-running completed work
- Is parameterised via `conf/analysis.config`; no paths are hardcoded in `main.nf`

### 2. All outputs redirected to this directory
The original scripts wrote outputs into `jingling_analysis/`. All outputs now land under `results/` within this directory. The only reads from outside this directory are the four external inputs listed below.

### 3. Updated scripts
Three scripts were updated or rewritten relative to Jingling's originals:

| Script | Change |
|---|---|
| `scripts/prepare_input_data_for_tgfm_updated.py` | Updated version with fixed NaN handling and allele alignment |
| `scripts/run_tgfm_updated.py` | Updated TGFM runner; imports `tgfm_nan_update_v2.py` instead of original `tgfm.py` |
| `scripts/tgfm_nan_update_v2.py` | Modified TGFM core with line 23 (`np.clip`) removed to fix a numerical instability |
| `scripts/make_gwas_summary_file.py` | Rewritten with `argparse`; discovers sumstat files dynamically rather than using hardcoded paths |
| `scripts/generate_window_per_chunk.py` | Rewritten with `argparse`; sorts windows genome-wide before chunking |

### 4. Robust handling of centromeric gaps
The original LD generation scripts would crash (exit code 13 from plink2) when a 3 Mb window fell entirely within a centromeric gap and contained no variants. The pipeline now detects exit code 13 and silently skips those windows rather than failing the entire chromosome job.

### 5. Fix for Nextflow config path-resolution bug
Nextflow evaluates `${params.xxx}` references *inside* a `params {}` block against the empty top-level map, producing literal `[:]` strings instead of the intended paths. All derived paths are now defined as Groovy variables *outside* the `params {}` block and referenced within it.

### 6. Smart strand-flip step
The pipeline checks at startup whether Jingling's pre-computed corrected-flip plink files already exist at `geno_corrected_flip_source` (set in `conf/analysis.config`). If they do, a lightweight `LINK_CORRECTED_FLIP` process symlinks them into `results/genotype/` and `CORRECT_FLIP` is skipped entirely. If they do not exist, `CORRECT_FLIP` runs the local `scripts/correct_flip.sh`, which uses `INPUT_PLINK_PREFIX` / `OUTPUT_PLINK_PREFIX` environment variables (no hardcoded paths) and delegates to the original Jingling script via `CORRECT_FLIP_DELEGATE`.

`INTERSECT_GWAS_GENO` runs after the flip step completes, reading from `params.plink_prefix` (the canonical corrected-flip prefix) in both cases.

### 7. Fixed `generate_windows.py` argument order
The script expects positional arguments `<bim_file> <output_file> <chromosome>`. The previous call passed four arguments in the wrong order; this is now corrected.

---

## External inputs (read-only)

| Input | Path |
|---|---|
| Human core exome genotype (raw) | `jingling_analysis/human_exome_core_genotype/humancoreexome_allchr_subset_included_in_ibd_analysis` |
| Pre-computed corrected-flip genotype | `jingling_analysis/human_core_exome_genotype/humancoreexome_filtered_plink_corrected_flip` |
| eQTL imputed genotype | `core_analysis_output/.../IBDverse_genotypes/2024_07_11-genotype_plate12345/plink` |
| eQTL summary statistics | `jingling_analysis/IBDverse_data/reformatted_eqtl_sumstat/` |
| GWAS summary statistics | `jingling_analysis/IBDverse_data/DeLange_sumstat/` |

Paths are set in `conf/analysis.config`.

---

## Running the pipeline

```bash
bsub < submit_nextflow_BH.sh
```

The submit script loads Singularity and Nextflow, then runs:

```bash
nextflow run workflow/main.nf \
    -profile sanger \
    -resume \
    -c nextflow.config \
    -c conf/analysis.config
```

To skip steps that have already completed (e.g. re-run only TGFM after LD matrices are done), set the corresponding skip flags in `conf/analysis.config`:

```groovy
skip_intersect     = true   // skip genotype filtering and GWAS intersect
skip_ld_generation = true   // skip window/LD/npy generation
skip_eqtl_fm       = true   // skip SuSiE eQTL fine-mapping
```

---

## Pipeline stages

See `summ.md` for a detailed description of each step. Briefly:

1. **Intersect & flip** — filter genotype and GWAS sumstats to a common SNP set; correct strand flips
2. **LD matrices** — generate per-window LD matrices (chr 1–22 in parallel), fix NaN entries, convert to `.npy`
3. **eQTL fine-mapping** — SuSiE fine-mapping per cell type × chromosome using `susie_eqtl_fine_mapping_for_tgfm.py`
4. **Summary files** — build tissue and GWAS lookup tables
5. **TGFM input** — integrate LD + eQTL posteriors + GWAS per chromosome; split into 50 genome-wide chunks
6. **TGFM** — fast run (no sampling) to estimate tissue priors; full run (100 samples) for final PIPs

Final results are in `results/tgfm_final/{CD,UC,IBD}/`.
