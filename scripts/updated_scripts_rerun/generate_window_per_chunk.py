import os
import glob
import numpy as np


def read_and_sort_tgfm_input_summaries(input_dir):
    # Find all summary files
    files = sorted(glob.glob(os.path.join(input_dir, "tgfm_chr*_input_data_summary.txt")))
    all_rows = []
    header = None

    for file in files:
        with open(file, 'r') as f:
            lines = f.readlines()
            if not lines:
                continue
            if header is None:
                header = lines[0].strip()
            for line in lines[1:]:
                line = line.strip()
                if line:
                    all_rows.append(line)

    # Sort by chromosome and start position extracted from window_name (first column)
    def sort_key(row):
        window_name = row.split('\t')[0]
        chrom, pos_range = window_name.split(':')
        start = int(pos_range.split('-')[0])
        return (int(chrom), start)

    all_rows = sorted(all_rows, key=sort_key)
    return header, all_rows


def split_into_chunks(data_lines, n_chunks):
    return np.array_split(data_lines, n_chunks)

def write_chunks_to_files(header, chunks, output_prefix):
    for i, chunk in enumerate(chunks):
        out_filename = f"{output_prefix}_chunk{i+1}_input_data_summary.txt"
        with open(out_filename, 'w') as out_file:
            if header:
                out_file.write(header + '\n')
            for line in chunk:
                out_file.write(line + '\n')

def main():
    input_dir = "/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/output_data/input_data_for_TGFM"  
    output_prefix = "tgfm"  
    n_chunks = 50

    header, all_rows = read_and_sort_tgfm_input_summaries(input_dir)
    print(f"Loaded {len(all_rows)} rows from TGFM files")
    chunks = split_into_chunks(all_rows, n_chunks)
    write_chunks_to_files(header, chunks, output_prefix)
    print(f"✅ Done! Split into {n_chunks} chunks with prefix '{output_prefix}_chunk#_input_data_summary.txt'.")

if __name__ == "__main__":
    main()
