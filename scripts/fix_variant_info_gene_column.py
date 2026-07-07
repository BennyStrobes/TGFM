"""
Fix existing combined variant-info TSV files that were written without the gene (ENSG) name
in column 0.  The original writer looped  `for gene, data in accumulated.items()`  but only
wrote  `\t'.join(row)`  (the 6-column bim row), producing chromosome number in col 0 instead
of the ENSG ID.

Reconstruction method:
  1. Read the gene_summary file to get the ordered list of Pass genes (same insertion order
     as `accumulated.items()` when the vi file was written).
  2. For each passing gene, load the npz and read `npz[gene + '__alpha'].shape[1]` to learn
     how many variant rows belong to that gene.
  3. Read the vi file sequentially; for each gene, consume N rows and prepend the gene name.

The fixed file is written in-place (old content replaced with corrected content).
"""

import os
import sys
import numpy as np


def fix_vi_file(vi_file, summary_file, npz_file):
    # --- ordered list of passing genes ---
    pass_genes = []
    with open(summary_file) as f:
        f.readline()  # skip header
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if parts[3] == 'Pass':
                pass_genes.append(parts[0])

    if not pass_genes:
        print(f"  SKIP (no passing genes): {vi_file}")
        return

    # --- variant counts per gene from npz ---
    npz = np.load(npz_file)
    gene_n_variants = {}
    for gene in pass_genes:
        key = gene + '__alpha'
        if key not in npz:
            print(f"  WARNING: {key} not in npz — skipping gene")
            gene_n_variants[gene] = 0
        else:
            gene_n_variants[gene] = npz[key].shape[1]
    npz.close()

    total_expected = sum(gene_n_variants.values())

    # --- read all data rows (skip header) ---
    with open(vi_file) as f:
        header = f.readline().rstrip('\n')
        data_rows = [line.rstrip('\n') for line in f]

    if len(data_rows) != total_expected:
        print(f"  ERROR: expected {total_expected} rows, found {len(data_rows)} — "
              f"skipping {vi_file}")
        return

    # Verify first column is a chromosome number (not already an ENSG ID)
    if data_rows and data_rows[0].startswith('ENSG'):
        print(f"  ALREADY FIXED: {vi_file}")
        return

    # --- rewrite with gene prefix ---
    fixed_rows = []
    idx = 0
    for gene in pass_genes:
        n = gene_n_variants[gene]
        for row in data_rows[idx: idx + n]:
            fixed_rows.append(gene + '\t' + row)
        idx += n

    with open(vi_file, 'w') as f:
        f.write(header + '\n')
        for row in fixed_rows:
            f.write(row + '\n')

    print(f"  FIXED ({len(fixed_rows)} rows, {len(pass_genes)} genes): {os.path.basename(vi_file)}")


def main(eqtl_dir):
    vi_files = sorted(f for f in os.listdir(eqtl_dir) if f.endswith('_variant_info.tsv'))
    print(f"Found {len(vi_files)} vi files to process")

    n_ok = n_skip = n_err = 0
    for vf in vi_files:
        stem = vf.replace('_variant_info.tsv', '')
        vi_path      = os.path.join(eqtl_dir, vf)
        summary_path = os.path.join(eqtl_dir, stem + '_gene_summary.txt')
        npz_path     = os.path.join(eqtl_dir, stem + '_susie.npz')

        if not os.path.isfile(summary_path):
            print(f"  SKIP (no summary): {vf}")
            n_skip += 1
            continue
        if not os.path.isfile(npz_path):
            print(f"  SKIP (no npz): {vf}")
            n_skip += 1
            continue

        try:
            fix_vi_file(vi_path, summary_path, npz_path)
            n_ok += 1
        except Exception as e:
            print(f"  ERROR processing {vf}: {e}")
            n_err += 1

    print(f"\nDone: {n_ok} fixed, {n_skip} skipped, {n_err} errors")


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python fix_variant_info_gene_column.py <eqtl_fine_mapping_dir>")
        sys.exit(1)
    main(sys.argv[1])
