import os
import re
import gzip
import csv
import glob
import argparse


def main():
    parser = argparse.ArgumentParser(description="Generate a GWAS summary file from intersected summary statistics.")
    parser.add_argument("--sumstats-dir", required=True,
                        help="Directory containing *geno_intersected_GWAS_sumstat.txt.gz files")
    parser.add_argument("--output-file", required=True,
                        help="Path to the output GWAS summary file")
    args = parser.parse_args()

    sample_size_column_index = 5

    trait_info = {}
    sample_size_dict = {}

    # get all the matching gwas sumstat file
    gwas_sumstat_pattern = os.path.join(args.sumstats_dir, "*geno_intersected_GWAS_sumstat.txt.gz")
    matching_file = glob.glob(gwas_sumstat_pattern)

    for file_path in matching_file:
        # Extract trait name from the GWAS summary statistics files
        filename = os.path.basename(file_path)
        trait_name = filename.split('_')[0]
        trait_info[trait_name] = file_path

        # Extract the sample size from GWAS summary statistics files
        try:
            with gzip.open(file_path, 'rt') as f:
                reader = csv.reader(f, delimiter='\t')
                header = next(reader, None)
                first_data_row = next(reader, None)
                if first_data_row and len(first_data_row) > sample_size_column_index:
                    sample_size_value = first_data_row[sample_size_column_index].strip()
                    sample_size_dict[trait_name] = sample_size_value
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            sample_size_dict[trait_name] = "NA"

    # Write the GWAS summary file
    with open(args.output_file, 'w') as out:
        out.write("trait_name\tsample_size\tsummary_statistics\n")
        for trait_name in sorted(trait_info.keys()):
            file_path = trait_info[trait_name]
            sample_size_value = sample_size_dict.get(trait_name, "NA")
            out.write(f"{trait_name}\t{sample_size_value}\t{file_path}\n")

    print(f"GWAS summary file written to: {args.output_file}")


if __name__ == "__main__":
    main()
