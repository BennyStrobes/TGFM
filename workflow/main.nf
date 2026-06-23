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
// Step 1a — Generate per-window PLINK beds and LD matrices (chr-parallel)
//
// Writes to params.plink_perwindow_dir and params.ld_matrix_dir.
// These directories are outside Nextflow's work dir; a sentinel signals
// completion so downstream processes can chain correctly.
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
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Step 1b — Re-compute LD for windows that contained NaN values
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
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Step 1c — Convert LD matrices to .npy and extract variant info
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
    echo "Finished chromosome ${chr}"
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Step 2 — Run TGFM per trait × chunk
// ─────────────────────────────────────────────────────────────────────────────
process RUN_TGFM {
    label 'mem_30GB_long'
    tag "${trait}_chunk${chunk}"
    container params.container_tgfm

    // TGFM writes results directly to params.output_dir/${trait};
    // no publishDir needed — the path is already the final destination.

    input:
    tuple val(trait), val(chunk)

    output:
    tuple val(trait), val(chunk)

    script:
    """
    set -euo pipefail

    mkdir -p "${params.output_dir}/${trait}"

    # TGFM expects to be run from Jingling's scripts directory
    cd "${params.tgfm_workdir}"

    echo "Running TGFM: trait=${trait}  chunk=${chunk}"

    python "${params.tgfm_script}" \\
        --trait-name "${trait}" \\
        --tgfm-input-data "${params.tgfm_input_dir}/tgfm_chunk${chunk}_input_data_summary.txt" \\
        --parallel-job-identifier "${chunk}" \\
        --tissue-summary-file "${params.tissue_summary_file}" \\
        --tissue-specific-prior "${params.tissue_prior_dir}/${trait}" \\
        --out "${params.output_dir}/${trait}"

    echo "Finished: trait=${trait}  chunk=${chunk}"
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Workflow
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    log.info """
    ╔══════════════════════════════════════════════════════╗
    ║           TGFM Pipeline                              ║
    ╠══════════════════════════════════════════════════════╣
    ║  Traits    : ${params.traits}
    ║  Chunks    : ${params.n_chunks}
    ║  Chromosomes: ${params.chromosomes}
    ║  Output    : ${params.output_dir}
    ╚══════════════════════════════════════════════════════╝
    """.stripIndent()

    // ── Step 1: LD matrix generation (parallelised per chromosome) ─────────
    if (!params.skip_ld_generation) {
        def chroms = parseChroms(params.chromosomes)
        chr_ch = Channel.from(chroms)

        // Chain the three sub-steps sequentially within each chromosome,
        // but all chromosomes run in parallel across the cluster.
        ld_done_ch   = GENERATE_LD_MATRIX(chr_ch)
        nan_done_ch  = FIX_NAN_LD_MATRIX(ld_done_ch)
        npy_done_ch  = GENERATE_NPY_VARIANT_INFO(nan_done_ch)

        // Gate TGFM on ALL chromosomes completing before any chunk starts
        all_chrs_done = npy_done_ch.collect()
    } else {
        // If LD generation is skipped, create a dummy completion signal
        all_chrs_done = Channel.value(['skipped'])
    }

    // ── Step 2: TGFM — fan out over traits × chunks ───────────────────────
    def traits     = params.traits.split(',').collect { it.trim() }
    def chunk_list = (1..params.n_chunks).toList()

    // Build trait-chunk pairs only after all chromosomes are done
    trait_chunk_ch = all_chrs_done
        .flatMap { _ ->
            def pairs = []
            traits.each { trait ->
                chunk_list.each { chunk -> pairs << [trait, chunk] }
            }
            pairs
        }

    RUN_TGFM(trait_chunk_ch)
}
