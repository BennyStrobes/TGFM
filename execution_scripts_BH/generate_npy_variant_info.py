import os
import numpy as np
import pandas as pd
import sys

def convert_ld_to_npy(ld_matrix_dir, npy_output_dir, window_name):
    """Convert LD matrix file in .vcor2 format (long lines with no terminators) to .npy"""
    vcor2_file = os.path.join(ld_matrix_dir, f"{window_name}.unphased.vcor2")
    vars_file = os.path.join(ld_matrix_dir, f"{window_name}.unphased.vcor2.vars")
    npy_file = os.path.join(npy_output_dir, f"{window_name}_ld.npy")

    if not os.path.exists(vcor2_file) or not os.path.exists(vars_file):
        print(f"Missing required files: {vcor2_file} or {vars_file}")
        return

    # Read the .vars file to determine the number of variants (N)
    with open(vars_file, "r") as f:
        variant_ids = f.read().splitlines()
    n_variants = len(variant_ids)

    # Read the .vcor2 file (long lines with no line terminators)
    with open(vcor2_file, "r") as f:
        # Read all lines (ignoring line terminators) and split into whitespace-separated elements
        elements = f.read().split()  # Split by any whitespace (spaces, tabs, newlines)
        
    # Convert to float32 array
    data = np.array(elements, dtype=np.float32)

    # Determine matrix shape
    expected_triangular = n_variants * (n_variants + 1) // 2
    expected_square = n_variants ** 2

    if len(data) == expected_triangular:
        # Case 1: Lower-triangular matrix
        matrix = np.zeros((n_variants, n_variants), dtype=np.float32)
        rows, cols = np.tril_indices(n_variants)
        matrix[rows, cols] = data
        matrix = matrix + matrix.T - np.diag(matrix.diagonal())  # Symmetrize
    elif len(data) == expected_square:
        # Case 2: Square matrix
        matrix = data.reshape(n_variants, n_variants)
    else:
        raise ValueError(
            f"Data size {len(data)} does not match triangular ({expected_triangular}) "
            f"or square ({expected_square}) matrix for {n_variants} variants."
        )

    # Save as .npy
    np.save(npy_file, matrix)
    print(f"Saved {npy_file}")



def convert_bim_to_variant_info(plink_file_dir, variant_info_output_dir, window_name):
    """Convert PLINK .bim file to a tab-delimited variant info file"""
    bim_file = os.path.join(plink_file_dir, f"{window_name}.bim")
    variant_info_file = os.path.join(variant_info_output_dir, f"{window_name}_variant_info.txt")

    if not os.path.exists(bim_file):
        print(f"BIM file not found: {bim_file}")
        return

    # Read BIM file (tab-delimited by default)
    df = pd.read_csv(bim_file, sep="\t", header=None)

    # Save as tab-delimited text file
    df.to_csv(variant_info_file, sep="\t", index=False, header=False)
    print(f"Saved {variant_info_file}")


if __name__ == "__main__":
    # Parse command-line arguments
    ld_matrix_dir = sys.argv[1]
    plink_file_dir = sys.argv[2]
    npy_output_dir = sys.argv[3]
    variant_info_output_dir = sys.argv[4]
    window_name = sys.argv[5]

    convert_ld_to_npy(ld_matrix_dir, npy_output_dir, window_name)
    convert_bim_to_variant_info(plink_file_dir, variant_info_output_dir, window_name)
