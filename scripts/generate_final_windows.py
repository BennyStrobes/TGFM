"""
Generate final window file for a chromosome by joining the basic window TSV with
the paths to the corresponding LD .npy and variant info files.

Usage:
    python generate_final_windows.py <chromosome> <window_input_dir> <npy_dir> <variant_info_dir> <output_dir>

Output:
    <output_dir>/window_file_chr<chromosome>.tsv
    Columns: window_name, chr, window_start, window_end, LD_matrix, variant_info_file
"""

import os
import sys
import pandas as pd


def process_chromosome(chromosome, window_input_dir, npy_dir, variant_info_dir, output_dir):
    window_input_file = os.path.join(window_input_dir, f"windows_chr{chromosome}.tsv")
    final_output      = os.path.join(output_dir,       f"window_file_chr{chromosome}.tsv")

    if not os.path.exists(window_input_file):
        print(f"Error: {window_input_file} not found.")
        sys.exit(1)

    df = pd.read_csv(window_input_file, sep="\t")
    df["LD_matrix"]         = df["window_name"].apply(lambda x: os.path.join(npy_dir,         f"{x}_ld.npy"))
    df["variant_info_file"] = df["window_name"].apply(lambda x: os.path.join(variant_info_dir, f"{x}_variant_info.txt"))

    os.makedirs(output_dir, exist_ok=True)
    df.to_csv(final_output, sep="\t", index=False)
    print(f"Chromosome {chromosome}: {len(df)} windows written to {final_output}")


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python generate_final_windows.py "
              "<chromosome> <window_input_dir> <npy_dir> <variant_info_dir> <output_dir>")
        sys.exit(1)

    process_chromosome(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
