"""
Scan all LD matrix windows for a given chromosome and write out the variant IDs
that appear in NaN-containing rows/columns.

Usage:
    python make_nan_indices.py <chromosome> <ld_matrix_dir> <output_file>

Arguments:
    chromosome      Integer chromosome number (e.g. 1)
    ld_matrix_dir   Directory containing <window>.unphased.vcor2[.vars] files
    output_file     TSV output: window_name <TAB> variant_id for every NaN variant
"""

import os
import sys
import glob
import numpy as np


def nan_variants_for_window(vcor2_file, vars_file):
    """Return list of variant IDs whose row/column in the LD matrix contains a NaN."""
    with open(vars_file) as f:
        variant_ids = f.read().splitlines()
    n = len(variant_ids)

    with open(vcor2_file) as f:
        elements = f.read().split()
    data = np.array(elements, dtype=np.float32)

    expected_tri = n * (n + 1) // 2
    expected_sq  = n * n

    if len(data) == expected_tri:
        matrix = np.zeros((n, n), dtype=np.float32)
        rows, cols = np.tril_indices(n)
        matrix[rows, cols] = data
        matrix = matrix + matrix.T - np.diag(matrix.diagonal())
    elif len(data) == expected_sq:
        matrix = data.reshape(n, n)
    else:
        print(f"WARNING: unexpected data size {len(data)} for {vcor2_file} "
              f"(expected triangular={expected_tri} or square={expected_sq}) — skipping")
        return []

    if not np.isnan(matrix).any():
        return []

    # Variants whose entire row or column is NaN
    all_nan_mask = np.isnan(matrix).all(axis=1) | np.isnan(matrix).all(axis=0)
    all_nan_idx  = np.where(all_nan_mask)[0]

    # Among remaining variants, those with any NaN
    keep_mask = ~all_nan_mask
    any_nan_idx = np.array([], dtype=int)
    if keep_mask.any():
        sub = matrix[keep_mask][:, keep_mask]
        any_nan_sub = np.isnan(sub).any(axis=1) | np.isnan(sub).any(axis=0)
        orig_idx    = np.where(keep_mask)[0]
        any_nan_idx = orig_idx[any_nan_sub]

    combined = np.unique(np.concatenate([all_nan_idx, any_nan_idx])).astype(int)
    print(f"  {os.path.basename(vcor2_file)}: {len(combined)} NaN variant(s)")
    return [variant_ids[i] for i in combined]


def main():
    if len(sys.argv) != 4:
        print("Usage: python make_nan_indices.py <chromosome> <ld_matrix_dir> <output_file>")
        sys.exit(1)

    chrom        = sys.argv[1]
    ld_matrix_dir = sys.argv[2]
    output_file  = sys.argv[3]

    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)

    pattern = os.path.join(ld_matrix_dir, f"{chrom}:*.unphased.vcor2")
    vcor2_files = sorted(glob.glob(pattern))

    if not vcor2_files:
        print(f"No .vcor2 files found for chromosome {chrom} in {ld_matrix_dir}")
        open(output_file, "w").close()
        return

    total = 0
    with open(output_file, "w") as out:
        for vcor2_file in vcor2_files:
            vars_file = vcor2_file + ".vars"
            if not os.path.exists(vars_file):
                print(f"WARNING: missing .vars file for {vcor2_file} — skipping")
                continue

            window_name = os.path.basename(vcor2_file).replace(".unphased.vcor2", "")
            nan_vars = nan_variants_for_window(vcor2_file, vars_file)
            for var_id in nan_vars:
                out.write(f"{window_name}\t{var_id}\n")
                total += 1

    print(f"Chromosome {chrom}: {total} NaN variant entries written to {output_file}")


if __name__ == "__main__":
    main()
