import pandas as pd
import os
import sys


def generate_windows(bim_file, output_file, chrom, window_size=3000000, step_size=1000000):
    """Generate overlapping genomic windows for a specific chromosome from a PLINK bim file."""
    if not os.path.exists(bim_file):
        print(f"Error: {bim_file} not found.")
        sys.exit(1)

    df = pd.read_csv(bim_file, sep=r'\s+', header=None, usecols=[0, 3], names=["chr", "pos"])
    df = df[df["chr"] == chrom]
    if df.empty:
        print(f"Warning: No data found for chromosome {chrom}.")
        return

    start_pos = max(1, df["pos"].min() - 1000000)
    last_variant_pos = df["pos"].max()

    windows = []
    pos = start_pos
    while pos <= last_variant_pos:
        window_start = pos
        window_end = pos + window_size
        windows.append([f"{chrom}:{window_start}-{window_end}", chrom, window_start, window_end])
        if last_variant_pos <= window_end:
            break
        pos += step_size

    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    window_df = pd.DataFrame(windows, columns=["window_name", "chr", "window_start", "window_end"])
    window_df.to_csv(output_file, sep="\t", index=False)
    print(f"Generated {len(windows)} windows: {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python generate_windows.py <bim_file> <output_file> <chromosome>")
        sys.exit(1)

    generate_windows(sys.argv[1], sys.argv[2], int(sys.argv[3]))
