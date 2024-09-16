import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os 
import pdb
import pickle
import argparse


def load_in_gwas_data_for_single_trait_first_trait(sumstat_file, chrom_num):
	g = open(sumstat_file)
	dictionary = {}
	variant_names = []
	variant_positions = []
	allele_1s = []
	allele_2s = []
	betas = []
	beta_vars = []

	head_count = 0
	for line in g:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Throw out variants not on this chromosome
		line_chrom = data[0]
		if line_chrom != chrom_num:
			continue

		# Extract relevent fields
		variant_name = data[1]
		variant_pos = int(data[2])
		a1 = data[3]
		a2 = data[4]
		beta = float(data[6])
		beta_var = float(data[7])

		# Add to array
		variant_names.append(variant_name)
		variant_positions.append(variant_pos)
		allele_1s.append(a1)
		allele_2s.append(a2)
		betas.append(beta)
		beta_vars.append(beta_var)
	g.close()

	dictionary['beta'] = np.asarray(betas)
	dictionary['beta_var'] = np.asarray(beta_vars)

	return dictionary, np.asarray(variant_names), np.asarray(variant_positions), np.asarray(allele_1s), np.asarray(allele_2s)

def load_in_gwas_data_for_single_trait_not_first_trait(sumstat_file, chrom_num, variant_names, variant_positions, allele_1, allele_2):
	g = open(sumstat_file)
	dictionary = {}
	betas = []
	beta_vars = []

	head_count = 0
	snp_counter = 0
	for line in g:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Throw out variants not on this chromosome
		line_chrom = data[0]
		if line_chrom != chrom_num:
			continue

		# Extract relevent fields
		variant_name = data[1]
		variant_pos = int(data[2])
		a1 = data[3]
		a2 = data[4]
		beta = float(data[6])
		beta_var = float(data[7])

		if variant_name != variant_names[snp_counter]:
			print('assumption error: GWAS summary statistic variants are not in the same order')
			pdb.set_trace()
		if variant_pos != variant_positions[snp_counter]:
			print('assumption error: GWAS summary statistic variants are not in the same order')
			pdb.set_trace()
		if a1 != allele_1[snp_counter]:
			print('assumption error: GWAS summary statistic variants are not in the same order')
			pdb.set_trace()
		if a2 != allele_2[snp_counter]:
			print('assumption error: GWAS summary statistic variants are not in the same order')
			pdb.set_trace()							

		# Add to array
		betas.append(beta)
		beta_vars.append(beta_var)

		snp_counter = snp_counter + 1
	g.close()

	dictionary['beta'] = np.asarray(betas)
	dictionary['beta_var'] = np.asarray(beta_vars)

	return dictionary


def load_in_gwas_data(gwas_summary_file, chrom_num):
	gwas_data_obj = {}
	gwas_data_obj['trait_sumstat'] = {}
	gwas_traits = []
	gwas_trait_sample_sizes = []
	f = open(gwas_summary_file)
	head_count = 0
	trait_count = 0

	# Loop through traits
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue

		# Parse data fields corresponding to this trait
		trait_name = data[0]
		trait_sample_size = int(data[1])
		sumstat_file = data[2]
		gwas_traits.append(trait_name)
		gwas_trait_sample_sizes.append(trait_sample_size)
		
		# Extract gwas summary statistics for this trait
		# Keep track of whether or not this is the first trait
		if trait_count == 0:
			single_trait_obj, variant_names, variant_positions, allele_1, allele_2 = load_in_gwas_data_for_single_trait_first_trait(sumstat_file, chrom_num)
		else:
			first_trait = False
			single_trait_obj = load_in_gwas_data_for_single_trait_not_first_trait(sumstat_file, chrom_num, variant_names, variant_positions, allele_1, allele_2)
		trait_count = trait_count + 1
		gwas_data_obj['trait_sumstat'][trait_name] = single_trait_obj
		
	f.close()

	# Add to global data object
	gwas_data_obj['gwas_trait_names'] = np.asarray(gwas_traits)
	gwas_data_obj['gwas_N'] = np.asarray(gwas_trait_sample_sizes)
	gwas_data_obj['variant_names'] = variant_names
	gwas_data_obj['variant_positions'] = variant_positions
	gwas_data_obj['variant_a1'] = allele_1
	gwas_data_obj['variant_a2'] = allele_2
	

	# Finally create dictionary mapping from variant name to variant index
	variant_name_to_variant_index_mapping = {}
	for ii, variant_name in enumerate(variant_names):
		variant_name_to_variant_index_mapping[variant_name] = ii
	gwas_data_obj['variant_to_index'] = variant_name_to_variant_index_mapping

	return gwas_data_obj





def load_in_gene_names_and_positions(tissue_summary_file, chrom_num):
	gene_data_obj = {}

	gene_tissue_names = []
	gene_names = []
	tissue_names = []
	gene_coords = []
	gene_variant_info_files = []
	susie_mu_files = []
	susie_mu_var_files = []
	susie_alpha_files = []
	susie_pmces_files = []

	head_count = 0
	f = open(tissue_summary_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')

		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue

		# Parse relevent fields corresponding to this line of the file
		tissue_name = data[0]
		tissue_eqtl_fine_mapping_output_file = data[2] + '_chr' + str(chrom_num) + '_gene_summary.txt'

		# Stream eqtl fine-mapping results file (for this tissue)
		g = open(tissue_eqtl_fine_mapping_output_file)
		head_count2 = 0
		for line2 in g:
			line2 = line2.rstrip()
			data2 = line2.split('\t')
			if head_count2 == 0:
				head_count2 = head_count2 + 1
				continue

			# Parse line
			gene_name = data2[0]
			line_chrom_num = data2[1]
			gene_coord = int(data2[2])
			gene_status = data2[3]

			# Skip genes that did not pass filters
			if gene_status != 'Pass':
				continue

			# Further parse line
			gene_variant_info_file = data2[4]
			gene_susie_alpha_file = data2[5]
			gene_susie_mu_file = data2[6]
			gene_susie_mu_var_file = data2[7]
			gene_susie_pmces_file = data2[8]
			gene_tissue_name = gene_name + '_' + tissue_name

			# Add information to global arrays
			gene_tissue_names.append(gene_tissue_name)
			gene_names.append(gene_name)
			tissue_names.append(tissue_name)
			gene_coords.append(gene_coord)
			gene_variant_info_files.append(gene_variant_info_file)
			susie_mu_files.append(gene_susie_mu_file)
			susie_mu_var_files.append(gene_susie_mu_var_file)
			susie_alpha_files.append(gene_susie_alpha_file)
			susie_pmces_files.append(gene_susie_pmces_file)
		g.close()
	f.close()


	# Add to gene data object
	gene_data_obj['gene_tissue_pairs'] = np.asarray(gene_tissue_names)
	gene_data_obj['genes'] = np.asarray(gene_names)
	gene_data_obj['tissues'] = np.asarray(tissue_names)
	gene_data_obj['gene_coord'] = np.asarray(gene_coords)
	gene_data_obj['gene_variant_info_files'] = np.asarray(gene_variant_info_files)
	gene_data_obj['susie_mu_files'] = np.asarray(susie_mu_files)
	gene_data_obj['susie_mu_var_files'] = np.asarray(susie_mu_var_files)
	gene_data_obj['susie_alpha_files'] = np.asarray(susie_alpha_files)
	gene_data_obj['susie_pmces_files'] = np.asarray(susie_pmces_files)

	return gene_data_obj




def load_in_variant_data(ld_variant_info_file):
	# initialize arrays to store data
	window_variants = []
	variant_positions = []
	variant_a1s = []
	variant_a2s = []


	# Stream variant info file
	g = open(ld_variant_info_file)
	for line2 in g:
		line2 = line2.rstrip()
		data2 = line2.split('\t')

		# parse line corresponding to a variant
		variant_name = data2[1]
		variant_pos = int(data2[3])
		variant_a1 = data2[4]
		variant_a2 = data2[5]

		# Add to arrays
		window_variants.append(variant_name)
		variant_positions.append(variant_pos)
		variant_a1s.append(variant_a1)
		variant_a2s.append(variant_a2)

	# Close file handle
	g.close()

	# Create mapping from variant name to window index
	window_variants = np.asarray(window_variants)
	mapping = {}
	for ii, variant in enumerate(window_variants):
		mapping[variant] = ii

	# Create mapping from variant name to alleles
	variant_to_allele_mapping = {}
	for ii, variant in enumerate(window_variants):
		variant_to_allele_mapping[variant] = (variant_a1s[ii], variant_a2s[ii])


	return np.asarray(window_variants), np.asarray(variant_positions), np.asarray(variant_a1s), np.asarray(variant_a2s), mapping, variant_to_allele_mapping

def create_gene_index_to_window_index_mapping(variant_to_window_index, variant_to_allele_mapping, gene_variant_info_file):
	gene_to_window = []
	flips = []

	g = open(gene_variant_info_file)

	for line2 in g:
		line2 = line2.rstrip()
		data2 = line2.split('\t')

		eqtl_variant_name = data2[1]
		eqtl_variant_pos = int(data2[3])
		eqtl_variant_a1 = data2[4]
		eqtl_variant_a2 = data2[5]

		gene_to_window.append(variant_to_window_index[eqtl_variant_name])

		gwas_a1, gwas_a2 = variant_to_allele_mapping[eqtl_variant_name]

		if gwas_a1 == eqtl_variant_a1 and gwas_a2 == eqtl_variant_a2:
			flips.append(1.0)
		elif gwas_a1 == eqtl_variant_a2 and gwas_a2 == eqtl_variant_a1:
			flips.append(-1.0)
		else:
			print('Errror in variant mapping: alleles dont line up')
			pdb.set_trace()
	g.close()
	return np.asarray(gene_to_window), np.asarray(flips)

def create_window_level_susie_no_flips(gene_local_to_global_mapping, local_susie_data, num_global_variants):
	global_susie_data = np.zeros((local_susie_data.shape[0], num_global_variants))
	global_susie_data[:, gene_local_to_global_mapping] = local_susie_data
	return global_susie_data, local_susie_data

def create_window_level_susie_w_flips(gene_local_to_global_mapping, gene_flips, local_susie_data, num_global_variants):
	for var_index, flip_value in enumerate(gene_flips):
		if flip_value == -1.0:
			local_susie_data[:, var_index] = local_susie_data[:, var_index]*-1.0
	global_susie_data, local_susie_data = create_window_level_susie_no_flips(gene_local_to_global_mapping, local_susie_data, num_global_variants)
	return global_susie_data, local_susie_data

def sample_eqtl_effects_from_susie_posterior_distribution(gene_susie_mu, gene_susie_alpha, gene_susie_mu_var):
	n_components = gene_susie_mu.shape[0]
	n_snps = gene_susie_mu.shape[1]

	sampled_eqtl_effects = np.zeros(n_snps)

	for component_iter in range(n_components):
		# Randomly draw snp index for component
		random_snp_index = np.random.choice(np.arange(n_snps).astype(int), replace=False, p=gene_susie_alpha[component_iter,:])
		
		effect_size_mean = gene_susie_mu[component_iter,random_snp_index]
		effect_size_var = gene_susie_mu_var[component_iter, random_snp_index]

		random_effect_size = np.random.normal(loc=effect_size_mean, scale=np.sqrt(effect_size_var))

		sampled_eqtl_effects[random_snp_index] = sampled_eqtl_effects[random_snp_index] + random_effect_size
	return sampled_eqtl_effects

def prepare_tgfm_data_in_a_single_window(window_name, window_start, window_end, ld_variant_info_file, LD, gene_data_obj, num_posterior_samples, cis_window_size):
	#############################
	# Load in window variant data
	window_variants, window_variant_positions, window_variant_a1s, window_variant_a2s, variant_to_window_index, variant_to_allele_mapping = load_in_variant_data(ld_variant_info_file)
	# Get number of variants in window
	n_window_variants = len(window_variants)

	##############################
	# Get gene tissue pairs in this window
	# We only include gene-tissue pairs who's tss is larger than (window-start + cis window size) and who's tss is smaller than (window-end - cis window size)
	# This is done so gene-model doesn't go outside window boundry
	# And is fine as we only consider results of genes within center MB
	window_gene_lb = window_start + cis_window_size
	window_gene_ub = window_end - cis_window_size
	window_gene_tissue_indices = (gene_data_obj['gene_coord'] >= window_gene_lb) & (gene_data_obj['gene_coord'] <= window_gene_ub)
	# Ignore windows with no genes
	if np.sum(window_gene_tissue_indices) == 0:
		return {}, False
	# Extract gene-tissue pairs in window
	window_gene_tissue_pairs = gene_data_obj['gene_tissue_pairs'][window_gene_tissue_indices]
	window_gene_names = gene_data_obj['genes'][window_gene_tissue_indices]
	window_tissue_names = gene_data_obj['tissues'][window_gene_tissue_indices]
	window_gene_variant_info_files = gene_data_obj['gene_variant_info_files'][window_gene_tissue_indices]
	window_gene_susie_mu_files = gene_data_obj['susie_mu_files'][window_gene_tissue_indices]
	window_gene_susie_mu_var_files = gene_data_obj['susie_mu_var_files'][window_gene_tissue_indices]
	window_gene_susie_alpha_files = gene_data_obj['susie_alpha_files'][window_gene_tissue_indices]


	# Initialize arrays to keep track of eqtl-susie data
	pmces_arr = []
	gene_variances = []
	posterior_samples_arr = []


	# Loop through gene-tissue pairs
	for gene_tissue_index, gene_tissue_name in enumerate(window_gene_tissue_pairs):

		# Load in SuSiE eQTL gene model for this gene-tissue pair
		gene_tissue_susie_mu = np.load(window_gene_susie_mu_files[gene_tissue_index])
		gene_tissue_susie_mu_var = np.load(window_gene_susie_mu_var_files[gene_tissue_index])
		gene_tissue_susie_alpha = np.load(window_gene_susie_alpha_files[gene_tissue_index])

		# Align gene variants with window variants
		gene_variant_info_file = window_gene_variant_info_files[gene_tissue_index]
		gene_index_to_window_index_mapping, gene_sign_flips = create_gene_index_to_window_index_mapping(variant_to_window_index, variant_to_allele_mapping, gene_variant_info_file)
		
		# Create window-level susie alpha
		window_susie_alpha, gene_tissue_susie_alpha = create_window_level_susie_no_flips(gene_index_to_window_index_mapping, gene_tissue_susie_alpha, n_window_variants)
		# Create window-level susie mu
		window_susie_mu, gene_tissue_susie_mu = create_window_level_susie_w_flips(gene_index_to_window_index_mapping, gene_sign_flips, gene_tissue_susie_mu, n_window_variants)
		# Create window-level susie mu_sd
		window_susie_mu_sd, gene_tissue_susie_mu_sd = create_window_level_susie_w_flips(gene_index_to_window_index_mapping, gene_sign_flips, np.sqrt(gene_tissue_susie_mu_var), n_window_variants)

		# Create susie PMCES
		gene_tissue_susie_pmces = np.sum(gene_tissue_susie_mu*gene_tissue_susie_alpha,axis=0)
		window_susie_pmces = np.sum(window_susie_mu*window_susie_alpha,axis=0)

		# Construct LD matrix limited to snps in cis window of gene-tissue pair
		gene_tissue_LD = LD[gene_index_to_window_index_mapping,:][:, gene_index_to_window_index_mapping]

		# Calculate variance of genetically predicted gene expression
		gene_variance  = np.dot(np.dot(gene_tissue_susie_pmces, gene_tissue_LD), gene_tissue_susie_pmces)

		# Get standardized PMCES for gene-tissue pair
		standardized_window_susie_pmces = window_susie_pmces/np.sqrt(gene_variance)

		# Store data in window-level array
		gene_variances.append(gene_variance)
		pmces_arr.append(standardized_window_susie_pmces)


		######################
		# Draw samples from susie eqtl posterior distribution
		######################
		n_gene_tissue_snps = len(gene_tissue_susie_pmces)
		sampled_effects = np.zeros((n_gene_tissue_snps, num_posterior_samples))
		# For number of samples, sample
		for sample_iter in range(num_posterior_samples):
			susie_sampled_eqtl_effects = sample_eqtl_effects_from_susie_posterior_distribution(gene_tissue_susie_mu, gene_tissue_susie_alpha, np.square(gene_tissue_susie_mu_sd))
			sampled_effects[:, sample_iter] = susie_sampled_eqtl_effects
		# Standardize and store data		
		sampled_gene_variances = np.diag(np.dot(np.dot(np.transpose(sampled_effects), gene_tissue_LD), sampled_effects))
		sampled_standardized_effects = sampled_effects/np.sqrt(sampled_gene_variances)
		posterior_samples_arr.append((sampled_standardized_effects, gene_index_to_window_index_mapping))
	

	######################
	# The sampled effects are very, very sparse
	# For sake of memory, save the samples in a sparse format
	# The sparse format is as follows:
	######################	
	sparse_sampled_gene_eqtl_pmces = []
	for sample_iter in range(num_posterior_samples):
		eqtl_mat = []
		for gene_itera, sampled_eqtl_tuple in enumerate(posterior_samples_arr):
			eqtl_gene_window = sampled_eqtl_tuple[0][:, sample_iter]
			eqtl_indices = sampled_eqtl_tuple[1]
			for ii, eqtl_effect in enumerate(eqtl_gene_window):
				if eqtl_effect == 0.0:
					continue
				eqtl_index = eqtl_indices[ii]
				eqtl_mat.append(np.asarray([gene_itera, eqtl_index, eqtl_effect]))
		eqtl_mat = np.asarray(eqtl_mat)
		sparse_sampled_gene_eqtl_pmces.append(eqtl_mat)

	# Get middle indices
	window_middle_start = window_start + 1000000
	window_middle_end = window_end - 1000000
	middle_gene_indices = np.where((gene_data_obj['gene_coord'][window_gene_tissue_indices] >= window_middle_start) & (gene_data_obj['gene_coord'][window_gene_tissue_indices] < window_middle_end))[0]
	middle_variant_indices = np.where((window_variant_positions >= window_middle_start) & (window_variant_positions < window_middle_end))[0]

	# Save to global object
	window_dictionary = {'gene_tissue_pairs': window_gene_tissue_pairs, 'genes': window_gene_names, 'tissues': window_tissue_names, 'gene_coord':gene_data_obj['gene_coord'][window_gene_tissue_indices], 'variants': window_variants, 'variant_positions': window_variant_positions, 'gene_eqtl_pmces': np.asarray(pmces_arr), 'gene_variances': np.asarray(gene_variances), 'sparse_sampled_gene_eqtl_pmces':sparse_sampled_gene_eqtl_pmces, 'middle_gene_indices': middle_gene_indices, 'middle_variant_indices': middle_variant_indices}

	return window_dictionary, True


# Convert gwas summary statistics to *STANDARDIZED* effect sizes
# Following SuSiE code found in these two places:
########1. https://github.com/stephenslab/susieR/blob/master/R/susie_rss.R  (LINES 277-279)
########2. https://github.com/stephenslab/susieR/blob/master/R/susie_ss.R (LINES 148-156 AND 203-205)
def convert_to_standardized_summary_statistics(gwas_beta_raw, gwas_beta_se_raw, gwas_sample_size, R, sigma2=1.0):
	gwas_z_raw = gwas_beta_raw/gwas_beta_se_raw

	XtX = (gwas_sample_size-1)*R
	Xty = np.sqrt(gwas_sample_size-1)*gwas_z_raw
	var_y = 1

	dXtX = np.diag(XtX)
	csd = np.sqrt(dXtX/(gwas_sample_size-1))
	csd[csd == 0] = 1

	XtX = (np.transpose((1/csd) * XtX) / csd)
	Xty = Xty / csd

	dXtX2 = np.diag(XtX)

	beta_scaled = (1/dXtX2)*Xty
	beta_se_scaled = np.sqrt(sigma2/dXtX2)

	return beta_scaled, beta_se_scaled

def extract_window_gwas_beta_and_gwas_beta_se(ld_variant_info_file, gwas_data_obj):
	gwas_beta = []
	gwas_beta_var = []

	g = open(ld_variant_info_file)
	for line2 in g:
		line2 = line2.rstrip()
		data2 = line2.split('\t')

		variant_name = data2[1]
		variant_a1 = data2[4]
		variant_a2 = data2[5]

		gwas_index = gwas_data_obj['variant_to_index'][variant_name]
		tmp_beta_arr = []
		tmp_beta_var_arr = []
		for trait_name in gwas_data_obj['gwas_trait_names']:
			tmp_beta_arr.append(gwas_data_obj['trait_sumstat'][trait_name]['beta'][gwas_index])
			tmp_beta_var_arr.append(gwas_data_obj['trait_sumstat'][trait_name]['beta_var'][gwas_index])
		gwas_beta.append(np.asarray(tmp_beta_arr))
		gwas_beta_var.append(np.asarray(tmp_beta_var_arr))
	g.close()

	gwas_beta = np.asarray(gwas_beta)
	gwas_beta_se = np.sqrt(np.asarray(gwas_beta_var))

	return np.transpose(gwas_beta), np.transpose(gwas_beta_se)



def convert_to_standardized_summary_statistics_wrapper(gwas_beta, gwas_beta_se, study_sample_sizes, LD):
	scaled_gwas_betas = []
	scaled_gwas_beta_ses = []
	for trait_iter, study_sample_size in enumerate(study_sample_sizes):
		
		scaled_gwas_beta, scaled_gwas_beta_se = convert_to_standardized_summary_statistics(gwas_beta[trait_iter,:], gwas_beta_se[trait_iter,:], float(study_sample_size), LD)
		
		scaled_gwas_betas.append(scaled_gwas_beta)
		scaled_gwas_beta_ses.append(scaled_gwas_beta_se)

	scaled_gwas_betas = np.asarray(scaled_gwas_betas)
	scaled_gwas_beta_ses = np.asarray(scaled_gwas_beta_ses)

	return scaled_gwas_betas, scaled_gwas_beta_ses







###########################
# Parse command line args
###########################
parser = argparse.ArgumentParser()
parser.add_argument('--chrom', default='None', type=str,
                    help='Chromosome number (e.g. 1)')
parser.add_argument('--cis-window-size', default=500000, type=int,
                    help='Size of cis-window')
parser.add_argument('--min-variants-per-window', default=50, type=int,
                    help='Only consider windows with at least this number of variants')
parser.add_argument('--window-file', default='None', type=str,
                    help='File overlapping windows')
parser.add_argument('--tissue-summary-file', default='None', type=str,
                    help='File containing information on each tissue and info on where that tissues gene models are stored')
parser.add_argument('--random-seed', default=None, type=int,
                    help='Random seed')
parser.add_argument('--gwas-summary-file', default='None', type=str,
                    help='File containing GWAS studies to include (note all gwas studies must include identical variants and ordering of said variants)')
parser.add_argument('--standardize-gwas-summary-statistics', default=False, action='store_true',
                    help='Boolean on whether or not to transform gwas summary statistics to version where genotype is standardized to mean 0, variance 1')
parser.add_argument('--num-posterior-samples', default=100, type=int,
                    help='Number of posterior samples to randomly draw from each gene-tissue pair')
parser.add_argument('--out', default='None', type=str,
                    help='Output stem to print tgfm input data to')
args = parser.parse_args()



# Set random seed if requested
if args.random_seed is not None:
	np.random.seed(args.random_seed)


######################
# Load in GWAS Data
######################
print('Loading in GWAS summary statistic data')
gwas_data_obj = load_in_gwas_data(args.gwas_summary_file, args.chrom)

######################
# Load in Gene names and positions
######################
print('Load in gene names and positions')
gene_data_obj = load_in_gene_names_and_positions(args.tissue_summary_file, args.chrom)



# Open outputful summarizing TGFM input (one line for each window)
tgfm_input_data_summary_file = args.out + '_chr' + args.chrom + '_input_data_summary.txt'
t = open(tgfm_input_data_summary_file,'w')
# Write header
t.write('window_name\tLD_npy_file\tTGFM_input_pkl\tTGFM_trait_input_pkl\n')


######################
# Loop through windows (and subsequently prepare data in each window)
######################
# Open window summary file handle
f = open(args.window_file)
# Used to skip header line
head_count = 0

# Stream file
for line in f:
	line = line.rstrip()
	data = line.split('\t')

	# Skip header
	if head_count == 0:
		head_count = head_count + 1
		continue

	# Skip windows not on this chromosome
	line_chrom = data[1]
	if line_chrom != args.chrom:
		continue

	# Parse line (corresponds to a single window)
	window_name = data[0]
	window_start_pos = int(data[2])
	window_end_pos = int(data[3])
	window_ld_file = data[4]  # Corresponds to matrix of dimension n_snps X n_snps
	ld_variant_info_file = data[5]  # Plink bim file corresponding to snps making up the ld file (and their ordering)

	print(window_name)

	# Load in LD
	LD = np.load(window_ld_file)

	# Get number of snps in window
	n_window_snps = LD.shape[0]
	# Skip windows with fewer snps than min-variants-per-window
	if n_window_snps < args.min_variants_per_window:
		print('Window ' + window_name + ' skipped because too few snps in window')
		continue

	# Prepare trait-agnostic TGFM data in a single window
	tgfm_data_obj, pass_filter_boolean = prepare_tgfm_data_in_a_single_window(window_name, window_start_pos, window_end_pos, ld_variant_info_file, LD, gene_data_obj, args.num_posterior_samples, args.cis_window_size)

	# Skip window because window has no genes
	if pass_filter_boolean == False:
		print('Window ' + window_name + ' skipped because it contained 0 gene-tissue pairs')
		continue

	# Extract GWAS Beta and gwas beta se for variants in this window
	gwas_beta, gwas_beta_se = extract_window_gwas_beta_and_gwas_beta_se(ld_variant_info_file, gwas_data_obj)

	# Standardize gwas effects
	if args.standardize_gwas_summary_statistics:
		gwas_beta, gwas_beta_se = convert_to_standardized_summary_statistics_wrapper(gwas_beta, gwas_beta_se, gwas_data_obj['gwas_N'], LD)

	# Add trait data
	tgfm_trait_data = {}
	tgfm_trait_data['gwas_beta'] = gwas_beta
	tgfm_trait_data['gwas_beta_se'] = gwas_beta_se
	tgfm_trait_data['gwas_sample_size'] = gwas_data_obj['gwas_N']
	tgfm_trait_data['gwas_study_names'] = gwas_data_obj['gwas_trait_names']


	# Save TGFM input data object
	tgfm_pkl_file = args.out + '_' + window_name + '_tgfm_input_data.pkl'
	g = open(tgfm_pkl_file, "wb")
	pickle.dump(tgfm_data_obj, g)
	g.close()

	# Save TGFM-trait input data object
	tgfm_trait_pkl_file = args.out + '_' + window_name + '_tgfm_input_trait_data.pkl'
	g = open(tgfm_trait_pkl_file, "wb")
	pickle.dump(tgfm_trait_data, g)
	g.close()

	# Print locations to summary file
	t.write(window_name + '\t' + window_ld_file + '\t' + tgfm_pkl_file + '\t' + tgfm_trait_pkl_file + '\n')
	t.flush()

# Close window summary file handle
f.close()
t.close()


