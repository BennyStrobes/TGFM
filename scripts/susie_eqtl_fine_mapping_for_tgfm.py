import sys
import gzip
import pandas as pd
import numpy as np
import os
import pdb
from pandas_plink import read_plink1_bin
import rpy2
import rpy2.robjects.numpy2ri as numpy2ri
import rpy2.robjects as ro
numpy2ri.activate()
from rpy2.robjects.packages import importr
susieR_pkg = importr('susieR')
import argparse


def _open(path, mode='r'):
    if path.endswith('.gz'):
        return gzip.open(path, mode + 't')
    return open(path, mode)


def extract_gwas_variants(gwas_sumstat_file, chrom_num):
	f = _open(gwas_sumstat_file)
	dictionary = {}
	head_count = 0  # Skip header
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Skip variants not on this chromosome
		if data[0] != chrom_num:
			continue
		# Rsid of this line
		rsid = data[1]
		# error check to make sure rsid not used twice in gwas
		if rsid in dictionary:
			print('ASSUMPTION ERROR: rsid found twice in gwas sumstat file')
			pdb.set_trace()
		dictionary[rsid] = 1
	f.close()

	return dictionary

def mean_impute_and_standardize_genotype(G_obj_geno):
	# Fill in missing values
	G_obj_geno_stand = np.copy(G_obj_geno)
	ncol = G_obj_geno_stand.shape[1]
	for col_iter in range(ncol):
		nan_indices = np.isnan(G_obj_geno[:,col_iter])
		non_nan_mean = np.mean(G_obj_geno[nan_indices==False, col_iter])
		G_obj_geno_stand[nan_indices, col_iter] = non_nan_mean


	# Standardize genotype (Mean 0, variance 1)
	G_obj_geno_stand = (G_obj_geno_stand -np.mean(G_obj_geno_stand,axis=0))/np.std(G_obj_geno_stand,axis=0)

	return G_obj_geno_stand

def save_gene_variant_info(gene_variant_info_output_file, rsids, positions, a0s, a1s, chroms):
	# Save in plink format
	t2 = open(gene_variant_info_output_file,'w')

	# Loop through variants
	for var_iter, rsid in enumerate(rsids):
		t2.write(chroms[var_iter] + '\t' + rsid + '\t' + '0' + '\t' + str(positions[var_iter]) + '\t' + a0s[var_iter] + '\t' + a1s[var_iter] + '\n')
	t2.close()
	return


def load_in_eqtl_genotype_data(genotype_stem, chrom_num, gwas_variants, filter_strand_ambiguous):
	# Load in genotype data across chromosome for eQTL data set
	genotype_stem = genotype_stem + chrom_num
	G_obj = read_plink1_bin(genotype_stem + '.bed', genotype_stem + '.bim', genotype_stem + '.fam', verbose=False)
	G_obj_geno = G_obj.values # Numpy 2d array of dimension num samples X num snps
	G_obj_chrom = np.asarray(G_obj.chrom)
	G_obj_pos = np.asarray(G_obj.pos)
	# a0 is the effect allele
	# For case of plink package, a0 is the first column in the plink bim file
	G_obj_a0 = np.asarray(G_obj.a0)
	G_obj_a1 = np.asarray(G_obj.a1)
	# RSids
	G_obj_rsids = np.asarray(G_obj.snp)
	# Sample names
	G_obj_sample_names = np.asarray(G_obj.fid + '_' + G_obj.iid)
	# Snp ids
	G_obj_snp_ids = 'chr' + G_obj_chrom + '_' + (G_obj_pos.astype(str)) + '_' + G_obj_a0 + '_' + G_obj_a1

	# Mean impute and standardize genotype
	G_obj_geno_stand = mean_impute_and_standardize_genotype(G_obj_geno)

	########################
	# Filter variants to those in gwas and (optionally to strand ambigious variants)
	valid_variants = [] # Boolean vector to keep track of which variants pass filters

	# Loop through variants
	for var_iter, var_rsid in enumerate(G_obj_rsids):
		# Initialize boolean for variant to True
		booler = True

		# Identify those variants not in the gwas
		if var_rsid not in gwas_variants:
			booler = False

		# Filter out those variants not on this chromosome
		if G_obj_chrom[var_iter] != chrom_num:
			booler = False

		# Identify strand ambigious variants
		if filter_strand_ambiguous:
			var_a0 = G_obj_a0[var_iter]
			var_a1 = G_obj_a1[var_iter]
			if var_a0 == 'A' and var_a1 == 'T':
				booler = False
			if var_a0 == 'T' and var_a1 == 'A':
				booler = False
			if var_a0 == 'C' and var_a1 == 'G':
				booler = False
			if var_a0 == 'G' and var_a1 == 'C':
				booler = False

		valid_variants.append(booler)
	# Numpy array of valid variants
	valid_variants = np.asarray(valid_variants)

	print('Number of variants before variant filtering: ' + str(len(valid_variants)))
	print('Number of variants after variant filtering: ' + str(int(np.sum(valid_variants))))

	# Now filter output data structures to these variants

	return G_obj_geno_stand[:, valid_variants], G_obj_rsids[valid_variants], G_obj_a0[valid_variants], G_obj_a1[valid_variants], G_obj_pos[valid_variants], G_obj_chrom[valid_variants], G_obj_sample_names

def load_in_per_gene_eqtl_summary_statistics(eqtl_sumstat_file, rsid_to_genotype_position, chrom_num, cis_window_size):
	eqtl_sumstat_obj = {}
	head_count = 0
	f = _open(eqtl_sumstat_file)
	n_malformed = 0
	# Loop through eqtl summary stats
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Skip malformed lines (expected 10 columns: GENE SNP CHR GENE_COORD SNP_BP A1 A2 N BETA BETA_VAR)
		if len(data) < 10:
			n_malformed += 1
			if n_malformed <= 5:
				print(f'WARNING: skipping malformed line ({len(data)} columns, expected 10): {repr(line[:120])}')
			continue
		# Throw out variant-gene pairs not on correct chromosome
		line_chrom_num = data[2]
		if line_chrom_num != chrom_num:
			continue

		# Throw out snps filtered out
		line_rsid = data[1]
		if line_rsid not in rsid_to_genotype_position:
			continue

		# Extract other relevent info
		ensamble_id = data[0]
		gene_coord = int(data[3])
		snp_bp = int(data[4])
		a1 = data[5]
		a2 = data[6]

		# Throw out snps not in cis region
		distance = np.abs(snp_bp - gene_coord)
		if distance > cis_window_size:
			continue

		# Extract more relevent info from line
		beta = float(data[8])
		beta_var = float(data[9])
		sample_size = int(data[7])

		# Add snp gene pair to eqtl sumstat object
		if ensamble_id not in eqtl_sumstat_obj:
			eqtl_sumstat_obj[ensamble_id] = {}
			eqtl_sumstat_obj[ensamble_id]['beta'] = []
			eqtl_sumstat_obj[ensamble_id]['beta_var'] = []
			eqtl_sumstat_obj[ensamble_id]['rsid'] = []
			eqtl_sumstat_obj[ensamble_id]['genotype_index'] = []
			eqtl_sumstat_obj[ensamble_id]['sample_size'] = -1
			eqtl_sumstat_obj[ensamble_id]['gene_coord'] = 0
		eqtl_sumstat_obj[ensamble_id]['beta'].append(beta)
		eqtl_sumstat_obj[ensamble_id]['beta_var'].append(beta_var)
		eqtl_sumstat_obj[ensamble_id]['rsid'].append(line_rsid)
		eqtl_sumstat_obj[ensamble_id]['genotype_index'].append(rsid_to_genotype_position[line_rsid])
		eqtl_sumstat_obj[ensamble_id]['sample_size'] = sample_size
		eqtl_sumstat_obj[ensamble_id]['gene_coord'] = gene_coord
	f.close()
	if n_malformed > 0:
		print(f'WARNING: skipped {n_malformed} malformed lines in {eqtl_sumstat_file}')

	return eqtl_sumstat_obj

def run_susie_eqtl_fine_mapping_with_eqtl_summary_stats(gwas_variants, eqtl_sumstat_file, chrom_num, genotype_stem, cis_window_size, filter_strand_ambiguous, min_cis_snps_per_gene, output_stem):
	#############################
	# Load in eQTL Genotype data
	#############################
	print('Load in Genotype data')
	G_mat, G_rsids, G_a0, G_a1, G_pos, G_chrom, G_sample_names = load_in_eqtl_genotype_data(genotype_stem, chrom_num, gwas_variants, filter_strand_ambiguous)

	#############################
	# Create dictionary mapping from rsid to genotype position
	#############################
	rsid_to_genotype_position = {}
	for ii, rsid in enumerate(G_rsids):
		rsid_to_genotype_position[rsid] = ii

	#############################
	# Load in eQTL Summary Statistics
	#############################
	print('Load in eQTL summary statistics')
	per_gene_eqtl_summary_statistics = load_in_per_gene_eqtl_summary_statistics(eqtl_sumstat_file, rsid_to_genotype_position, chrom_num, cis_window_size)

	#############################
	# Fit gene model in each gene, independently; accumulate results in memory
	#############################
	print('Fit SuSiE eQTL gene models')
	summary_rows = []  # [gene, chrom, coord, status, combined_susie_or_NA, combined_vi_or_NA]
	accumulated = {}   # gene -> {alpha, mu, mu_var, pmces, variant_rows}

	genes = np.asarray([*per_gene_eqtl_summary_statistics])
	for gene in genes:

		# Extract relevent info for gene
		gene_beta = np.asarray(per_gene_eqtl_summary_statistics[gene]['beta'])
		gene_beta_var = np.asarray(per_gene_eqtl_summary_statistics[gene]['beta_var'])
		gene_rsids = np.asarray(per_gene_eqtl_summary_statistics[gene]['rsid'])
		gene_indices = np.asarray(per_gene_eqtl_summary_statistics[gene]['genotype_index'])
		eqtl_sample_size = per_gene_eqtl_summary_statistics[gene]['sample_size']
		gene_position = per_gene_eqtl_summary_statistics[gene]['gene_coord']

		# Get genotype data corresponding to this gene
		gene_geno = G_mat[:, gene_indices]
		gene_geno_rsids = G_rsids[gene_indices]
		gene_geno_a0 = G_a0[gene_indices]
		gene_geno_a1 = G_a1[gene_indices]
		gene_geno_snp_pos = G_pos[gene_indices]
		gene_chroms = G_chrom[gene_indices]

		#  Ignore genes with no or very few cis snps
		if len(gene_beta) < min_cis_snps_per_gene:
			print('gene skipped because it contained 0 or very small number of cis snps')
			summary_rows.append([gene, chrom_num, str(gene_position), 'Fail_too_few_snps', 'NA', 'NA'])
			continue

		# Error check to make sure rsids line up
		if np.array_equal(gene_geno_rsids, gene_rsids) == False:
			print('Assumption error: variants dont match between sumstats and genotype')
			pdb.set_trace()

		# Compute LD matrix across gene snps
		LD = np.corrcoef(np.transpose(gene_geno))

		# Run eQTL variant fine-mapping with SuSiE
		susie_fitted = susieR_pkg.susie_rss(bhat=gene_beta.reshape((len(gene_beta),1)), shat=(np.sqrt(gene_beta_var)).reshape((len(gene_beta_var),1)), R=LD, n=eqtl_sample_size, L=10)

		# Test whether are 0 identified susie components for this gene
		if type(susie_fitted.rx2('sets').rx2('cs_index')).__name__ == 'NULLType':
			summary_rows.append([gene, chrom_num, str(gene_position), 'Fail_purity_filter', 'NA', 'NA'])
			continue

		# This gene has passed purity filter — extract posteriors
		alpha = np.asarray(susie_fitted.rx2('alpha'))
		mu    = np.asarray(susie_fitted.rx2('mu'))
		mu2   = np.asarray(susie_fitted.rx2('mu2'))
		pmces = np.sum(alpha * mu, axis=0)
		mu_var = mu2 - np.square(mu)

		# Collect variant rows (one row per cis-SNP: chr, rsid, cm, pos, a0, a1)
		variant_rows = [
			[gene_chroms[i], gene_rsids[i], '0', str(int(gene_geno_snp_pos[i])), gene_geno_a0[i], gene_geno_a1[i]]
			for i in range(len(gene_rsids))
		]
		accumulated[gene] = {
			'alpha': alpha, 'mu': mu, 'mu_var': mu_var, 'pmces': pmces,
			'variant_rows': variant_rows,
		}
		summary_rows.append([gene, chrom_num, str(gene_position), 'Pass', None, None])

	# Combined output file paths (one npz + one tsv per chromosome × cell type)
	combined_npz_file = output_stem + '_chr' + str(chrom_num) + '_susie.npz'
	combined_vi_file  = output_stem + '_chr' + str(chrom_num) + '_variant_info.tsv'

	if accumulated:
		# Save all genes' posteriors into one compressed npz
		npz_dict = {}
		for gene, data in accumulated.items():
			npz_dict[gene + '__alpha']  = data['alpha']
			npz_dict[gene + '__mu']     = data['mu']
			npz_dict[gene + '__mu_var'] = data['mu_var']
			npz_dict[gene + '__pmces']  = data['pmces']
		np.savez_compressed(combined_npz_file, **npz_dict)

		# Save all genes' variant info into one tsv with a leading gene column
		with open(combined_vi_file, 'w') as vi:
			vi.write('gene\tchr\trsid\tcm\tpos\ta0\ta1\n')
			for gene, data in accumulated.items():
				for row in data['variant_rows']:
					vi.write('\t'.join(row) + '\n')

	# Write gene summary (2 file columns instead of 5)
	output_summary_file = output_stem + '_chr' + str(chrom_num) + '_gene_summary.txt'
	with open(output_summary_file, 'w') as t:
		t.write('Gene\tCHR\tGENE_COORD\tINFO\tcombined_susie_file\tcombined_variant_info_file\n')
		for row in summary_rows:
			if row[3] == 'Pass':
				row[4] = combined_npz_file
				row[5] = combined_vi_file
			t.write('\t'.join(row) + '\n')

	return

def run_susie_eqtl_fine_mapping_with_individual_data(gwas_variants, expression_file, chrom_num, genotype_stem, cis_window_size, filter_strand_ambiguous, min_cis_snps_per_gene, output_stem):
	#############################
	# Load in eQTL Genotype data
	#############################
	print('Load in Genotype data')
	G_mat, G_rsids, G_a0, G_a1, G_pos, G_chrom, G_sample_names = load_in_eqtl_genotype_data(genotype_stem, chrom_num, gwas_variants, filter_strand_ambiguous)

	#############################
	# Fit gene model in each gene, independently; accumulate results in memory
	#############################
	print('Fit SuSiE eQTL gene models')
	summary_rows = []  # [gene, chrom, coord, status, combined_susie_or_NA, combined_vi_or_NA]
	accumulated = {}   # gene -> {alpha, mu, mu_var, pmces, variant_rows}

	f = _open(expression_file)
	head_count = 0  # To identify header
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Header
		if head_count == 0:
			head_count = head_count + 1
			# Make sure expression sample names and genotype sample names match
			expr_sample_names = np.asarray(data[3:])
			if np.array_equal(G_sample_names, expr_sample_names) == False:
				print('ASSUMPTION ERROR: Sample names in genotype file do not match sample names in expression matrix')
				pdb.set_trace()
			continue

		# Standard line (corresponding to a gene)
		ensamble_id = data[0]
		gene_chrom_num = data[1]
		# Skip genes not on this chromosome
		if gene_chrom_num != chrom_num:
			continue
		gene_position = int(data[2])
		expr_vec = np.asarray(data[3:]).astype(float)

		# Get indices of variants corresponding to cis window of this gene
		cis_window_start = gene_position - cis_window_size
		cis_window_end = gene_position + cis_window_size
		cis_snp_indices = (G_pos >= cis_window_start) & (G_pos < cis_window_end)

		#  Ignore genes with no or very few cis snps
		if np.sum(cis_snp_indices) < min_cis_snps_per_gene:
			print('gene skipped because it contained 0 or very small number of cis snps')
			summary_rows.append([ensamble_id, gene_chrom_num, str(gene_position), 'Fail_too_few_snps', 'NA', 'NA'])
			continue

		# Extract standardized matrix of cis snps around the gene
		gene_geno = G_mat[:, cis_snp_indices]

		# Run eQTL variant fine-mapping with SuSiE
		susie_fitted = susieR_pkg.susie(gene_geno, expr_vec, L=10)

		# Test whether are 0 identified susie components for this gene
		if type(susie_fitted.rx2('sets').rx2('cs_index')).__name__ == 'NULLType':
			summary_rows.append([ensamble_id, gene_chrom_num, str(gene_position), 'Fail_purity_filter', 'NA', 'NA'])
			continue

		# This gene has passed purity filter — extract posteriors
		alpha = np.asarray(susie_fitted.rx2('alpha'))
		mu    = np.asarray(susie_fitted.rx2('mu'))
		mu2   = np.asarray(susie_fitted.rx2('mu2'))
		pmces = np.sum(alpha * mu, axis=0)
		mu_var = mu2 - np.square(mu)

		# Collect variant rows (one row per cis-SNP: chr, rsid, cm, pos, a0, a1)
		cis_rsids  = G_rsids[cis_snp_indices]
		cis_pos    = G_pos[cis_snp_indices]
		cis_a0     = G_a0[cis_snp_indices]
		cis_a1     = G_a1[cis_snp_indices]
		cis_chroms = G_chrom[cis_snp_indices]
		variant_rows = [
			[cis_chroms[i], cis_rsids[i], '0', str(int(cis_pos[i])), cis_a0[i], cis_a1[i]]
			for i in range(len(cis_rsids))
		]
		accumulated[ensamble_id] = {
			'alpha': alpha, 'mu': mu, 'mu_var': mu_var, 'pmces': pmces,
			'variant_rows': variant_rows,
		}
		summary_rows.append([ensamble_id, gene_chrom_num, str(gene_position), 'Pass', None, None])

	f.close()

	# Combined output file paths (one npz + one tsv per chromosome × cell type)
	combined_npz_file = output_stem + '_chr' + str(chrom_num) + '_susie.npz'
	combined_vi_file  = output_stem + '_chr' + str(chrom_num) + '_variant_info.tsv'

	if accumulated:
		# Save all genes' posteriors into one compressed npz
		npz_dict = {}
		for gene, data in accumulated.items():
			npz_dict[gene + '__alpha']  = data['alpha']
			npz_dict[gene + '__mu']     = data['mu']
			npz_dict[gene + '__mu_var'] = data['mu_var']
			npz_dict[gene + '__pmces']  = data['pmces']
		np.savez_compressed(combined_npz_file, **npz_dict)

		# Save all genes' variant info into one tsv with a leading gene column
		with open(combined_vi_file, 'w') as vi:
			vi.write('gene\tchr\trsid\tcm\tpos\ta0\ta1\n')
			for gene, data in accumulated.items():
				for row in data['variant_rows']:
					vi.write('\t'.join(row) + '\n')

	# Write gene summary (2 file columns instead of 5)
	output_summary_file = output_stem + '_chr' + str(chrom_num) + '_gene_summary.txt'
	with open(output_summary_file, 'w') as t:
		t.write('Gene\tCHR\tGENE_COORD\tINFO\tcombined_susie_file\tcombined_variant_info_file\n')
		for row in summary_rows:
			if row[3] == 'Pass':
				row[4] = combined_npz_file
				row[5] = combined_vi_file
			t.write('\t'.join(row) + '\n')

	return






###########################
# Parse command line args
###########################
parser = argparse.ArgumentParser()
parser.add_argument('--eqtl-data-type', default='Individual', type=str,
                    help='Type of eqtl data used (either Individual or SumStat)')
parser.add_argument('--expression-file', default='None', type=str,
                    help='Expression matrix file name')
parser.add_argument('--chrom', default='None', type=str,
                    help='Chromosome number (e.g. 1)')
parser.add_argument('--genotype-stem', default='None', type=str,
                    help='Stem of genotype plink file')
parser.add_argument('--cis-window-size', default=500000, type=int,
                    help='Size of cis-window')
parser.add_argument('--min-cis-snps-per-gene', default=10, type=int,
                    help='Minimum number of cis snps per gene')
parser.add_argument('--filter-strand-ambiguous', default=False, action='store_true',
                    help='Filter out strand ambiguous variants')
parser.add_argument('--eqtl-sumstat', default='None', type=str,
                    help='File containing eQTL summary statistics.')
parser.add_argument('--gwas-sumstat', default='None', type=str,
                    help='File containing GWAS summary statistics. This is used as we only include variants in a gene model if those variants are found in our GWAS study')
parser.add_argument('--out', default='None', type=str,
                    help='Output stem to print eqtl fine-mapping results to')
args = parser.parse_args()


##########################
# First extract dictionary list of variants on this chromsome used in the GWAS
##########################
gwas_variants = extract_gwas_variants(args.gwas_sumstat, args.chrom)


##########################
# Run eQTL fine-mapping
##########################
if args.eqtl_data_type == 'Individual':  # Individual level data
	run_susie_eqtl_fine_mapping_with_individual_data(gwas_variants, args.expression_file, args.chrom, args.genotype_stem, args.cis_window_size, args.filter_strand_ambiguous, args.min_cis_snps_per_gene, args.out)
elif args.eqtl_data_type == 'SumStat':  # Summary statistics
	run_susie_eqtl_fine_mapping_with_eqtl_summary_stats(gwas_variants, args.eqtl_sumstat, args.chrom, args.genotype_stem, args.cis_window_size, args.filter_strand_ambiguous, args.min_cis_snps_per_gene, args.out)
else:
	print('ASSUMPTION ERROR: ' + str(args.eqtl_data_type) + ' is not a valid eqtl data type')
	pdb.set_trace()
