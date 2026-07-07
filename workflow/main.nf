#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
// Helper: parse a chromosome range/list string ("1-22", "1,5,10", "1-5,22")
// ─────────────────────────────────────────────────────────────────────────────
def parseChroms(String spec) {
    def chroms = [] as Set
    spec.split(',').each { part ->
        def m = part.trim() =~ /^(\d+)-(\d+)$/
        if (m) {
            (m[0][1].toInteger()..m[0][2].toInteger()).each { chroms << it }
        } else {
            chroms << part.trim().toInteger()
        }
    }
    return chroms.sort()
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 0 — Intersect GWAS sumstats with genotype (depends on corrected genotype)
//
// Creates:
//   results/genotype/humancoreexome_filtered_plink.{bed,bim,fam}
//   results/sumstats/{IBD,CD,UC}_geno_intersected_GWAS_sumstat.txt.gz
// ─────────────────────────────────────────────────────────────────────────────
process INTERSECT_GWAS_GENO {
    label 'mem_100GB'
    tag "intersect_gwas_geno"
    container params.container_tgfm

    input:
    val _ready  // depends on CORRECT_FLIP

    output:
    val true

    script:
    """
    set -euo pipefail

    GENO_PREFIX="${params.plink_prefix}"
    SUMSTAT_DIR="${params.gwas_sumstat_dir}"
    THREE_TRAIT_SNPS="${params.three_trait_snps}"
    GENO_OUT_DIR="${params.geno_output_dir}"
    SUMSTATS_OUT_DIR="${params.sumstats_output_dir}"

    mkdir -p "\${GENO_OUT_DIR}" "\${SUMSTATS_OUT_DIR}"

    WORK_DIR="\${GENO_OUT_DIR}/intersect_work"
    mkdir -p "\${WORK_DIR}"

    FILTERED_PLINK_PREFIX="\${GENO_OUT_DIR}/humancoreexome_filtered_plink"

    # ── Build the final SNP list: (IBD ∩ CD ∩ UC) ∩ Genotype ────────────────
    cut -f2 "\${GENO_PREFIX}.bim" | sort > "\${WORK_DIR}/geno_snps.txt"
    sort "\${THREE_TRAIT_SNPS}" > "\${WORK_DIR}/trait_snps.txt"

    comm -12 \
        "\${WORK_DIR}/trait_snps.txt" \
        "\${WORK_DIR}/geno_snps.txt" \
        > "\${WORK_DIR}/final_common_snps.txt"

    echo "Final SNP count:"
    wc -l "\${WORK_DIR}/final_common_snps.txt"

    # ── Filter genotype to common SNPs using plink2 ──────────────────────────
    "${params.plink2_bin}" \
        --bfile "\${GENO_PREFIX}" \
        --extract "\${WORK_DIR}/final_common_snps.txt" \
        --make-bed \
        --out "\${FILTERED_PLINK_PREFIX}"

    # ── Filter all three GWAS sumstats to common SNPs ────────────────────────
    for trait in IBD CD UC; do
        zcat "\${SUMSTAT_DIR}/\${trait}_intersected_GWAS_sumstat.txt.gz" | \
        awk -v snpfile="\${WORK_DIR}/final_common_snps.txt" '
            BEGIN {
                while ((getline line < snpfile) > 0) keep[line] = 1
            }
            NR==1 || (\$2 in keep) { print }
        ' | gzip -c > "\${SUMSTATS_OUT_DIR}/\${trait}_geno_intersected_GWAS_sumstat.txt.gz"
        echo "Filtered sumstat written: \${SUMSTATS_OUT_DIR}/\${trait}_geno_intersected_GWAS_sumstat.txt.gz"
    done

    echo "INTERSECT_GWAS_GENO done."
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 1a — Symlink pre-computed corrected-flip files into results/genotype/
//
// Used instead of CORRECT_FLIP when geno_corrected_flip_source exists.
// Creates symlinks at params.plink_prefix so all downstream steps can find them.
// ─────────────────────────────────────────────────────────────────────────────
process LINK_CORRECTED_FLIP {
    label 'mem_low'
    tag "link_corrected_flip"

    output:
    val true

    script:
    """
    set -euo pipefail

    mkdir -p "${params.geno_output_dir}"

    for ext in bed bim fam; do
        src="${params.geno_corrected_flip_source}.\${ext}"
        dst="${params.plink_prefix}.\${ext}"
        if [[ ! -e "\${dst}" ]]; then
            ln -sf "\${src}" "\${dst}"
            echo "Linked \${dst} -> \${src}"
        else
            echo "Already exists: \${dst}"
        fi
    done

    echo "LINK_CORRECTED_FLIP done."
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 1b — Correct strand flips from scratch (runs when source does not exist)
//
// INPUT_PLINK_PREFIX / OUTPUT_PLINK_PREFIX / CORRECT_FLIP_DELEGATE are passed
// as env vars so scripts/correct_flip.sh contains no hardcoded paths.
// ─────────────────────────────────────────────────────────────────────────────
process CORRECT_FLIP {
    label 'mem_100GB'
    tag "correct_flip"
    // No container — the flip script may call system plink1

    output:
    val true

    script:
    """
    set -euo pipefail

    export INPUT_PLINK_PREFIX="${params.geno_raw_prefix}"
    export OUTPUT_PLINK_PREFIX="${params.plink_prefix}"
    export PLINK2_BIN="${params.plink2_bin}"
    export FLIP_ALLELES_FILE="${params.flip_alleles_file}"

    mkdir -p "\$(dirname "\${OUTPUT_PLINK_PREFIX}")"

    # Skip if output already exists (idempotent on -resume)
    if [[ -f "\${OUTPUT_PLINK_PREFIX}.bed" && -f "\${OUTPUT_PLINK_PREFIX}.bim" && -f "\${OUTPUT_PLINK_PREFIX}.fam" ]]; then
        echo "Corrected flip files already exist — skipping"
        ls -lh "\${OUTPUT_PLINK_PREFIX}".{bed,bim,fam}
        echo "CORRECT_FLIP done (skipped)."
        exit 0
    fi

    bash "${params.correct_flip_script}"

    ls -lh "\${OUTPUT_PLINK_PREFIX}.bed"
    echo "CORRECT_FLIP done."
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 2 — Generate window files per chromosome
// ─────────────────────────────────────────────────────────────────────────────
process GENERATE_WINDOWS {
    label 'mem_100GB'
    tag "windows_chr${chr}"
    container params.container_python

    input:
    val chr

    output:
    val chr

    script:
    """
    set -euo pipefail

    mkdir -p "${params.window_output_dir}"

    python "${params.generate_windows_script}" \
        "${params.plink_prefix}.bim" \
        "${params.window_output_dir}/windows_chr${chr}.tsv" \
        ${chr}

    echo "GENERATE_WINDOWS chr${chr} done."
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 3 — Generate per-window PLINK beds and LD matrices (chr-parallel)
//
// Writes to params.plink_perwindow_dir and params.ld_matrix_dir.
// PLINK exit code 13 = no variants in window → skip (e.g. centromeric gap).
// ─────────────────────────────────────────────────────────────────────────────
process GENERATE_LD_MATRIX {
    label 'mem_100GB'
    tag "ld_chr${chr}"

    input:
    val chr

    output:
    val chr

    script:
    """
    set -euo pipefail

    plink2_bin="${params.plink2_bin}"
    #chmod +x "\${plink2_bin}"
    "\${plink2_bin}" --version

    window_file="${params.window_dir}/windows_chr${chr}.tsv"
    plink_out="${params.plink_perwindow_dir}"
    ld_out="${params.ld_matrix_dir}"

    mkdir -p "\${plink_out}" "\${ld_out}"

    while IFS=\$'\\t' read -r window_name chr_col start end; do
        [[ "\${window_name}" == "window_name" ]] && continue

        # Extract per-window PLINK bed
        # Exit code 13 = no variants remaining (e.g. centromeric gap) — skip window
        make_bed_rc=0
        "\${plink2_bin}" --bfile "${params.plink_prefix}" \\
            --chr "\${chr_col}" \\
            --from-bp "\${start}" \\
            --to-bp "\${end}" \\
            --make-bed \\
            --out "\${plink_out}/\${window_name}" || make_bed_rc=\$?
        if [[ \${make_bed_rc} -eq 13 ]]; then
            echo "No variants in window \${window_name} — skipping"
            continue
        elif [[ \${make_bed_rc} -ne 0 ]]; then
            exit \${make_bed_rc}
        fi

        # Compute LD matrix
        "\${plink2_bin}" --bfile "\${plink_out}/\${window_name}" \\
            --r2-unphased square 'ref-based' \\
            --out "\${ld_out}/\${window_name}"

    done < "\${window_file}"

    echo "GENERATE_LD_MATRIX chr${chr} done."
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 4a — Identify NaN-containing windows per chromosome
// ─────────────────────────────────────────────────────────────────────────────
process MAKE_NAN_INDICES {
    label 'mem_100GB'
    tag "nan_idx_chr${chr}"
    container params.container_python

    input:
    val chr

    output:
    val chr

    script:
    """
    set -euo pipefail

    mkdir -p "${params.nan_indices_dir}"

    python "${params.make_nan_indices_script}" \
        ${chr} \
        "${params.ld_matrix_dir}" \
        "${params.nan_indices_dir}/nan_indices_chr${chr}.txt"

    echo "MAKE_NAN_INDICES chr${chr} done."
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 4b — Concatenate per-chromosome NaN index files into a single file
// ─────────────────────────────────────────────────────────────────────────────
process COLLECT_NAN_INDICES {
    label 'mem_low'
    tag "collect_nan"
    container params.container_python

    input:
    val all_chrs   // list of all chr values (collected)

    output:
    val true

    script:
    """
    set -euo pipefail

    NAN_DIR="${params.nan_indices_dir}"
    OUT_FILE="${params.nan_indices_file}"

    # Concatenate all per-chr nan index files
    cat "\${NAN_DIR}"/nan_indices_chr*.txt > "\${OUT_FILE}" 2>/dev/null || true

    echo "Total NaN entries:"
    wc -l "\${OUT_FILE}" || echo "0"
    echo "COLLECT_NAN_INDICES done."
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 5 — Re-compute LD for windows that contained NaN values
// ─────────────────────────────────────────────────────────────────────────────
process FIX_NAN_LD_MATRIX {
    label 'mem_100GB'
    tag "nan_chr${chr}"

    input:
    val chr

    output:
    val chr

    script:
    """
    set -euo pipefail

    plink2_bin="${params.plink2_bin}"
    exclude_file="${params.nan_indices_file}"
    window_file="${params.window_dir}/windows_chr${chr}.tsv"
    plink_out="${params.plink_perwindow_dir}"
    ld_out="${params.ld_matrix_dir}"

    # Build lookup: window_name -> space-separated variant list to exclude
    declare -A WINDOW_WITH_NAN
    while IFS=\$'\\t' read -r window variant; do
        WINDOW_WITH_NAN["\${window}"]="\${WINDOW_WITH_NAN["\${window}"]:-} \${variant}"
    done < "\${exclude_file}"

    while IFS=\$'\\t' read -r window_name chr_col start end; do
        [[ "\${window_name}" == "window_name" ]] && continue

        # Skip windows without NaN variants
        [[ -z "\${WINDOW_WITH_NAN["\${window_name}"]:-}" ]] && continue

        echo "Fixing NaN window: \${window_name}"

        exclude_tmp=\$(mktemp)
        for v in \${WINDOW_WITH_NAN["\${window_name}"]}; do
            echo "\${v}"
        done > "\${exclude_tmp}"

        make_bed_rc=0
        "\${plink2_bin}" --bfile "${params.plink_prefix}" \\
            --chr "\${chr_col}" \\
            --from-bp "\${start}" \\
            --to-bp "\${end}" \\
            --exclude "\${exclude_tmp}" \\
            --make-bed \\
            --out "\${plink_out}/\${window_name}" || make_bed_rc=\$?
        rm "\${exclude_tmp}"
        if [[ \${make_bed_rc} -eq 13 ]]; then
            echo "No variants in window \${window_name} after exclusions — skipping"
            continue
        elif [[ \${make_bed_rc} -ne 0 ]]; then
            exit \${make_bed_rc}
        fi

        "\${plink2_bin}" --bfile "\${plink_out}/\${window_name}" \\
            --r2-unphased square 'ref-based' \\
            --out "\${ld_out}/\${window_name}"

    done < "\${window_file}"

    echo "FIX_NAN_LD_MATRIX chr${chr} done."
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 6 — Convert LD matrices to .npy and extract variant info
// ─────────────────────────────────────────────────────────────────────────────
process GENERATE_NPY_VARIANT_INFO {
    label 'mem_100GB'
    tag "npy_chr${chr}"

    container params.container_python

    input:
    val chr

    output:
    val chr

    script:
    """
    set -euo pipefail

    ld_matrix_dir="${params.ld_matrix_dir}"
    plink_dir="${params.plink_perwindow_dir}"
    npy_dir="${params.npy_output_dir}"
    variant_info_dir="${params.variant_info_output_dir}"

    mkdir -p "\${npy_dir}" "\${variant_info_dir}"

    echo "Processing all windows for chromosome ${chr}..."
    found=0

    for vcor2_file in "\${ld_matrix_dir}/${chr}":*.unphased.vcor2; do
        [[ -f "\${vcor2_file}" ]] || continue
        found=1
        window_name=\$(basename "\${vcor2_file}" .unphased.vcor2)
        echo "  Processing window: \${window_name}"
        python "${params.generate_npy_script}" \\
            "\${ld_matrix_dir}" \\
            "\${plink_dir}" \\
            "\${npy_dir}" \\
            "\${variant_info_dir}" \\
            "\${window_name}"
    done

    [[ "\${found}" -eq 0 ]] && echo "WARNING: no .vcor2 files found for chr ${chr}"
    echo "GENERATE_NPY_VARIANT_INFO chr${chr} done."
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 7 — Generate final window files merging paths to LD/npy/variant info
// ─────────────────────────────────────────────────────────────────────────────
process GENERATE_FINAL_WINDOWS {
    label 'mem_100GB'
    tag "finalwin_chr${chr}"
    container params.container_python

    input:
    val chr

    output:
    val chr

    script:
    """
    set -euo pipefail

    python "${params.generate_final_windows_script}" \
        ${chr} \
        "${params.window_output_dir}" \
        "${params.npy_output_dir}" \
        "${params.variant_info_output_dir}" \
        "${params.window_output_dir}"

    echo "GENERATE_FINAL_WINDOWS chr${chr} done."
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 7b — Discover available cell types from the eQTL sumstat directory.
//
// Runs on a compute node (NFS-mounted) and emits one line per cell type:
//   file_prefix <TAB> folder_path
// e.g. "B_ct\t/path/to/dMean__B_ct_all"
// ─────────────────────────────────────────────────────────────────────────────
process DISCOVER_CELL_TYPES {
    label 'mem_low'
    tag "discover_cell_types"
    container params.container_tgfm

    output:
    path "cell_types.txt"

    script:
    """
    set -euo pipefail

    EQTL_SUMSTAT_BASE="${params.eqtl_sumstat_dir}"
    OUT="cell_types.txt"
    > "\${OUT}"

    for FOLDER in "\${EQTL_SUMSTAT_BASE}"/dMean__*_ct_all; do
        [[ -d "\${FOLDER}" ]] || continue
        folder_name=\$(basename "\${FOLDER}")

        # Strip leading dMean__ and trailing _ct_all to get the middle part
        # e.g. dMean__B_ct_all -> B   or  dMean__T_Tregs_ct_all -> T_Tregs
        inner=\${folder_name#dMean__}
        inner=\${inner%_ct_all}

        FILE_PREFIX="\${inner}_ct"
        echo "\${FILE_PREFIX}\t\${FOLDER}" >> "\${OUT}"
    done

    if [[ ! -s "\${OUT}" ]]; then
        echo "ERROR: No dMean__*_ct_all folders found in \${EQTL_SUMSTAT_BASE}" >&2
        exit 1
    fi

    echo "Discovered cell types:"
    cat "\${OUT}"
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 8 — SuSiE eQTL fine-mapping, one job per (chromosome × cell type)
//
// Output: results/eqtl_fine_mapping/{file_prefix}_chr{chr}_gene_summary.txt
//         + associated npy files
// ─────────────────────────────────────────────────────────────────────────────
process SUSIE_EQTL_FM {
    label 'mem_100GB_long'
    tag "susie_chr${chr}_${file_prefix}"
    container params.container_tgfm

    input:
    tuple val(chr), val(file_prefix), val(eqtl_folder)

    output:
    tuple val(chr), val(file_prefix)

    script:
    """
    set -euo pipefail

    EQTL_SUMSTAT="${eqtl_folder}/reformatted_cis_nominal1_eqtl.${chr}.tsv"
    EQTL_OUTPUT_STEM="${params.eqtl_fm_dir}/${file_prefix}"
    PLINK_GENO_FILE_STEM="${params.eqtl_plink_dir}/imputed_chr"
    GWAS_SUMSTAT="${params.sumstats_output_dir}/IBD_geno_intersected_GWAS_sumstat.txt.gz"

    mkdir -p "${params.eqtl_fm_dir}"

    echo "Processing: chr${chr}, cell type: ${file_prefix}"

    python "${params.susie_eqtl_script}" \\
        --eqtl-data-type SumStat \\
        --chrom "${chr}" \\
        --genotype-stem "\${PLINK_GENO_FILE_STEM}" \\
        --eqtl-sumstat "\${EQTL_SUMSTAT}" \\
        --gwas-sumstat "\${GWAS_SUMSTAT}" \\
        --filter-strand-ambiguous \\
        --out "\${EQTL_OUTPUT_STEM}"

    echo "SUSIE_EQTL_FM chr${chr} ${file_prefix} done."
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 9 — Build tissue summary file from eQTL fine-mapping outputs
//
// Scans results/eqtl_fine_mapping/ for *_ct_chr1_gene_summary.txt files,
// extracts cell type names, and writes results/tissue_summary_file.txt.
// ─────────────────────────────────────────────────────────────────────────────
process MAKE_TISSUE_SUMMARY {
    label 'mem_low'
    tag "tissue_summary"
    container params.container_python

    input:
    val _ready   // collected signal from all SUSIE_EQTL_FM jobs

    output:
    val true

    script:
    """
    set -euo pipefail

    python3 - <<'PYEOF'
import os
import glob

eqtl_fm_dir = "${params.eqtl_fm_dir}"
output_file = "${params.tissue_summary_file}"

# Find all chr1 gene summary files to identify available tissues
pattern = os.path.join(eqtl_fm_dir, "*_ct_chr1_gene_summary.txt")
files = sorted(glob.glob(pattern))

if not files:
    raise RuntimeError(f"No *_ct_chr1_gene_summary.txt files found in {eqtl_fm_dir}")

os.makedirs(os.path.dirname(output_file), exist_ok=True)

rows = []
for f in files:
    basename = os.path.basename(f)
    # Remove _ct_chr1_gene_summary.txt suffix to get cell_type prefix
    # e.g. "T_0_ct_chr1_gene_summary.txt" -> "T_0"
    suffix = "_ct_chr1_gene_summary.txt"
    if basename.endswith(suffix):
        cell_type = basename[: -len(suffix)]
    else:
        # Fallback: strip known suffix pattern
        cell_type = basename.replace("_chr1_gene_summary.txt", "")

    # output_stem = eqtl_fm_dir/{cell_type} (without _ct suffix per TGFM convention)
    # The _ct suffix is part of the file naming convention; tissue_summary uses stem without _ct
    if cell_type.endswith("_ct"):
        stem_name = cell_type[:-3]  # strip trailing _ct
    else:
        stem_name = cell_type
    output_stem = os.path.join(eqtl_fm_dir, stem_name)

    # Get sample size from the gene summary file (first data row, assume column index 1 or NA)
    sample_size = "NA"
    try:
        with open(f) as fh:
            header = fh.readline()
            first = fh.readline().strip()
            if first:
                parts = first.split("\\t")
                # sample_size column varies; use NA as safe default
                sample_size = "NA"
    except Exception:
        pass

    rows.append((cell_type, sample_size, output_stem))

with open(output_file, "w") as out:
    out.write("tissue_name\\tsample_size\\toutput_stem\\n")
    for tissue_name, sample_size, output_stem in rows:
        out.write(f"{tissue_name}\\t{sample_size}\\t{output_stem}\\n")

print(f"Tissue summary written: {output_file} ({len(rows)} tissues)")
PYEOF

    echo "MAKE_TISSUE_SUMMARY done."
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 10 — Build GWAS summary file from intersected sumstats
// ─────────────────────────────────────────────────────────────────────────────
process MAKE_GWAS_SUMMARY {
    label 'mem_low'
    tag "gwas_summary"
    container params.container_python

    input:
    val _ready   // signal from INTERSECT_GWAS_GENO

    output:
    val true

    script:
    """
    set -euo pipefail

    mkdir -p \$(dirname "${params.gwas_summary_file}")

    python "${params.make_gwas_summary_script}" \
        --sumstats-dir "${params.sumstats_output_dir}" \
        --output-file "${params.gwas_summary_file}"

    echo "MAKE_GWAS_SUMMARY done."
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 11 — Prepare TGFM input data per chromosome
// ─────────────────────────────────────────────────────────────────────────────
process PREPARE_INPUT {
    label 'mem_100GB'
    tag "prepare_input_chr${chr}"
    container params.container_tgfm

    input:
    val chr

    output:
    val chr

    script:
    """
    set -euo pipefail

    mkdir -p "${params.tgfm_input_raw_dir}"

    python "${params.prepare_input_script}" \
        --chrom ${chr} \
        --window-file "${params.window_dir}/window_file_chr${chr}.tsv" \
        --tissue-summary-file "${params.tissue_summary_file}" \
        --gwas-summary-file "${params.gwas_summary_file}" \
        --standardize-gwas-summary-statistics \
        --cis-window-size 1000000 \
        --out "${params.tgfm_input_raw_dir}/tgfm"

    echo "PREPARE_INPUT chr${chr} done."
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 12 — Split TGFM input into N chunks for parallel processing
// ─────────────────────────────────────────────────────────────────────────────
process GENERATE_WINDOW_PER_CHUNK {
    label 'mem_low'
    tag "chunk_split"
    container params.container_python

    input:
    val _ready   // collected signal from all PREPARE_INPUT jobs

    output:
    val true

    script:
    """
    set -euo pipefail

    mkdir -p "${params.tgfm_input_chunks_dir}"

    python "${params.generate_window_chunk_script}" \
        --input-dir "${params.tgfm_input_raw_dir}" \
        --output-dir "${params.tgfm_input_chunks_dir}" \
        --n-chunks ${params.n_chunks}

    echo "GENERATE_WINDOW_PER_CHUNK done."
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 13 — Run TGFM without sampling (50 chunks × 3 traits)
// ─────────────────────────────────────────────────────────────────────────────
process TGFM_WITHOUT_SAMPLING {
    label 'mem_100GB'
    tag "${trait}_chunk${chunk}"
    container params.container_tgfm

    input:
    tuple val(trait), val(chunk)

    output:
    tuple val(trait), val(chunk)

    script:
    """
    set -euo pipefail

    mkdir -p "${params.tgfm_no_sample_dir}/${trait}"

    python "${params.run_tgfm_no_sample_script}" \
        --trait-name "${trait}" \
        --tgfm-input-data "${params.tgfm_input_chunks_dir}/tgfm_chunk${chunk}_input_data_summary.txt" \
        --parallel-job-identifier "${chunk}" \
        --out "${params.tgfm_no_sample_dir}/${trait}" \
        --n-components 10 \
        --p-value-threshold 0.05 \
        --gene-tissue-pip-threshold 0.2

    echo "TGFM_WITHOUT_SAMPLING: trait=${trait}  chunk=${chunk} done."
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 14 — Compute tissue-specific prior per trait
// ─────────────────────────────────────────────────────────────────────────────
process TISSUE_PRIOR {
    label 'mem_100GB'
    tag "tissue_prior_${trait}"
    container params.container_tgfm

    input:
    // tuple: (trait, list-of-chunks-done)
    tuple val(trait), val(_chunks_done)

    output:
    val trait

    script:
    """
    set -euo pipefail

    mkdir -p "${params.tissue_prior_output_dir}/${trait}"

    # Write the chunk identifier file (1..n_chunks, one per line)
    CHUNK_ID_FILE=\$(mktemp)
    for i in \$(seq 1 ${params.n_chunks}); do
        echo "\${i}"
    done > "\${CHUNK_ID_FILE}"

    python "${params.tissue_prior_script}" \
        --trait-name "${trait}" \
        --tissue-summary-file "${params.tissue_summary_file}" \
        --tgfm-parallel-job-identifier-file "\${CHUNK_ID_FILE}" \
        --tgfm-without-sampling-output "${params.tgfm_no_sample_dir}/${trait}" \
        --out "${params.tissue_prior_output_dir}/${trait}"

    rm "\${CHUNK_ID_FILE}"
    echo "TISSUE_PRIOR trait=${trait} done."
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 15 — Final TGFM run (50 chunks × 3 traits)
// ─────────────────────────────────────────────────────────────────────────────
process FINAL_TGFM {
    label 'mem_30GB_long'
    tag "${trait}_chunk${chunk}"
    container params.container_tgfm

    input:
    tuple val(trait), val(chunk)

    output:
    tuple val(trait), val(chunk)

    script:
    """
    set -euo pipefail

    mkdir -p "${params.tgfm_final_dir}/${trait}"

    # run_tgfm_updated.py imports tgfm_nan_update_v2 from the same directory
    cd "${params.bh_dir}/scripts"

    python "${params.run_tgfm_updated_script}" \
        --trait-name "${trait}" \
        --tgfm-input-data "${params.tgfm_input_chunks_dir}/tgfm_chunk${chunk}_input_data_summary.txt" \
        --parallel-job-identifier "${chunk}" \
        --tissue-summary-file "${params.tissue_summary_file}" \
        --tissue-specific-prior "${params.tissue_prior_dir}/${trait}" \
        --out "${params.tgfm_final_dir}/${trait}"

    echo "FINAL_TGFM: trait=${trait}  chunk=${chunk} done."
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Workflow
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    log.info """
    ╔══════════════════════════════════════════════════════╗
    ║           TGFM Pipeline (full redesign)              ║
    ╠══════════════════════════════════════════════════════╣
    ║  Traits       : ${params.traits}
    ║  Chunks       : ${params.n_chunks}
    ║  Chromosomes  : ${params.chromosomes}
    ║  Output       : ${params.output_dir}
    ║  MCP mode     : ${params.mcp}
    ║  eQTL sumstat : ${params.eqtl_sumstat_dir}
    ║  skip_intersect    : ${params.skip_intersect}
    ║  skip_ld_generation: ${params.skip_ld_generation}
    ║  skip_eqtl_fm      : ${params.skip_eqtl_fm}
    ╚══════════════════════════════════════════════════════╝
    """.stripIndent()

    def chroms     = parseChroms(params.chromosomes)
    def traits     = params.traits.split(',').collect { it.trim() }
    def chunk_list = (1..params.n_chunks).toList()

    chr_ch = Channel.from(chroms)

    // ── STEP 1: Provide corrected-flip genotype ──────────────────────────────
    // If the pre-computed source files exist, symlink them (fast, no re-running).
    // Otherwise run CORRECT_FLIP from scratch.
    def _src = params.geno_corrected_flip_source
    if (_src && file("${_src}.bim").exists()) {
        log.info "Pre-computed corrected-flip found at ${_src} — symlinking, skipping CORRECT_FLIP"
        flip_done_ch = LINK_CORRECTED_FLIP()
    } else {
        log.info "No pre-computed corrected-flip found — running CORRECT_FLIP"
        flip_done_ch = CORRECT_FLIP()
    }

    // ── STEP 0: Intersect GWAS with corrected genotype (depends on flip step) ─
    if (!params.skip_intersect) {
        intersect_done_ch = INTERSECT_GWAS_GENO(flip_done_ch)
    } else {
        intersect_done_ch = Channel.value(true)
    }

    // ── STEPS 2-7: LD / window generation ────────────────────────────────────
    if (!params.skip_ld_generation) {
        // Step 2: generate window files (depends on correct_flip being done)
        windows_done_ch = GENERATE_WINDOWS(
            flip_done_ch.flatMap { chroms }
        )

        // Step 3: generate LD matrices (depends on window files)
        ld_done_ch = GENERATE_LD_MATRIX(windows_done_ch)

        // Step 4a: identify NaN windows per chromosome
        nan_idx_done_ch = MAKE_NAN_INDICES(ld_done_ch)

        // Step 4b: collect all NaN indices into single file
        all_nan_done_ch = COLLECT_NAN_INDICES(nan_idx_done_ch.collect())

        // Step 5: fix NaN LD matrices — needs all nan indices collected first
        // Broadcast the all_nan signal across all chromosomes
        nan_ready_chr_ch = all_nan_done_ch.flatMap { chroms }
        fix_done_ch      = FIX_NAN_LD_MATRIX(nan_ready_chr_ch)

        // Step 6: convert LD to npy
        npy_done_ch = GENERATE_NPY_VARIANT_INFO(fix_done_ch)

        // Step 7: generate final window files
        final_win_done_ch = GENERATE_FINAL_WINDOWS(npy_done_ch)

        // Gate downstream on all chromosomes completing
        all_ld_done = final_win_done_ch.collect()
    } else {
        all_ld_done = Channel.value(['skipped'])
    }

    // ── STEP 10: GWAS summary file (depends on intersect, runs in parallel with LD) ──
    gwas_summary_done_ch = MAKE_GWAS_SUMMARY(intersect_done_ch)

    // ── STEP 8: SuSiE eQTL fine-mapping (one job per chromosome × cell type) ───
    if (!params.skip_eqtl_fm) {
        // Discover cell types on a compute node (NFS-mounted), then split into
        // one (file_prefix, folder_path) pair per line of the output file.
        cell_types_ch = DISCOVER_CELL_TYPES()
            .splitText()
            .map { line ->
                def parts = line.trim().split('\t')
                [parts[0], parts[1]]   // [file_prefix, eqtl_folder]
            }

        // Combine: every chr with every cell type, gated on intersect completion
        susie_pairs_ch = intersect_done_ch
            .flatMap { chroms }
            .combine(cell_types_ch)
            .map { chr, file_prefix, eqtl_folder -> [chr, file_prefix, eqtl_folder] }

        susie_done_ch  = SUSIE_EQTL_FM(susie_pairs_ch)
        all_susie_done = susie_done_ch.collect()
    } else {
        all_susie_done = Channel.value(['skipped'])
    }

    // ── STEP 9: Tissue summary file ───────────────────────────────────────────
    tissue_summary_done_ch = MAKE_TISSUE_SUMMARY(all_susie_done)

    // ── STEP 11: Prepare TGFM input per chromosome ───────────────────────────
    // Depends on: LD done, tissue summary done, gwas summary done
    prepare_ready_ch = all_ld_done
        .combine(tissue_summary_done_ch)
        .combine(gwas_summary_done_ch)

    prepare_chr_ch    = prepare_ready_ch.flatMap { chroms }
    prepare_done_ch   = PREPARE_INPUT(prepare_chr_ch)
    all_prepare_done  = prepare_done_ch.collect()

    // ── STEP 12: Split into chunks ────────────────────────────────────────────
    chunk_split_done_ch = GENERATE_WINDOW_PER_CHUNK(all_prepare_done)

    // ── STEP 13: TGFM without sampling (50 chunks × 3 traits) ─────────────────
    no_sample_pairs_ch = chunk_split_done_ch.flatMap { _ ->
        def pairs = []
        traits.each { trait ->
            chunk_list.each { chunk -> pairs << [trait, chunk] }
        }
        pairs
    }
    no_sample_done_ch = TGFM_WITHOUT_SAMPLING(no_sample_pairs_ch)

    // ── STEP 14: Tissue prior — per trait, after all chunks for that trait done ─
    // Group completed (trait, chunk) pairs by trait, then run TISSUE_PRIOR per trait
    per_trait_done_ch = no_sample_done_ch
        .groupTuple(by: 0, size: chunk_list.size())
    tissue_prior_done_ch = TISSUE_PRIOR(per_trait_done_ch)

    // ── STEP 15: Final TGFM (50 chunks × 3 traits) ───────────────────────────
    // For each trait, fan back out to all chunks once the tissue prior is done
    final_pairs_ch = tissue_prior_done_ch.flatMap { trait ->
        chunk_list.collect { chunk -> [trait, chunk] }
    }
    FINAL_TGFM(final_pairs_ch)
}
