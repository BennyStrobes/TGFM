import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import numpy as np 
import os 
import pdb
import pickle
import tgfm
import rpy2
import rpy2.robjects.numpy2ri as numpy2ri
import rpy2.robjects as ro
import scipy.stats
ro.conversion.py2ri = numpy2ri
numpy2ri.activate()
from rpy2.robjects.packages import importr
susieR_pkg = importr('susieR')
import argparse


def extract_tissue_names(tissue_summary_file):
	# Extract ordered array of tissue names from tissue summary file
	f = open(tissue_summary_file)
	ordered_tissue_names = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		ordered_tissue_names.append(data[0])
	f.close()
	# Convert from python list to numpy array
	ordered_tissue_names = np.asarray(ordered_tissue_names)

	# Create dictionary mapping from tissue name to tissue position/index
	tissue_to_position_mapping = {}
	for i, val in enumerate(ordered_tissue_names):
		tissue_to_position_mapping[val] = i

	return ordered_tissue_names, tissue_to_position_mapping



def extract_full_gene_variant_ld(standardized_eqtl_effects, variant_ld):
	gene_variant_ld = np.dot(standardized_eqtl_effects,variant_ld) # Ngenes X n_variants
	expression_covariance = np.dot(gene_variant_ld, np.transpose(standardized_eqtl_effects))
	np.fill_diagonal(expression_covariance, 1.0)
	dd = np.diag(1.0/np.sqrt(np.diag(expression_covariance)))
	ge_ld = np.dot(np.dot(dd, expression_covariance),dd)
	top = np.hstack((ge_ld, gene_variant_ld))
	bottom = np.hstack((np.transpose(gene_variant_ld), variant_ld))
	full_ld = np.vstack((top,bottom))
	return full_ld




def extract_middle_genetic_elements(ordered_genes, middle_gene_indices, ordered_variants, middle_variant_indices):
	dicti = {}
	for gene_name in ordered_genes[middle_gene_indices]:
		dicti[gene_name] = 1
	for variant_name in ordered_variants[middle_variant_indices.astype(int)]:
		dicti[variant_name] = 1
	return dicti


def update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_obj):
	# Alphas
	susie_alpha = susie_obj.rx2('alpha')
	tgfm_obj.alpha_phi = susie_alpha[:,:(tgfm_obj.G)]
	tgfm_obj.beta_phi = susie_alpha[:,(tgfm_obj.G):]

	# Mus
	susie_mu = susie_obj.rx2('mu')
	tgfm_obj.alpha_mu = susie_mu[:,:(tgfm_obj.G)]
	tgfm_obj.beta_mu = susie_mu[:,(tgfm_obj.G):]

	# susie_mu_var
	susie_mu_var = susie_obj.rx2('mu2') - np.square(susie_mu)
	tgfm_obj.alpha_var = susie_mu_var[:,:(tgfm_obj.G)]
	tgfm_obj.beta_var = susie_mu_var[:,(tgfm_obj.G):]

	return tgfm_obj


def elbo_calc(z_vec, LD, samp_size, alpha, mu, mu2, KL_terms):
	bb = alpha*mu
	b_bar = np.sum(bb,axis=0)
	postb2 = alpha*mu2
	elbo_term1 = samp_size -1
	elbo_term2 = -2.0*np.sum(np.sqrt(samp_size-1)*b_bar*z_vec)
	elbo_term3 = np.sum(b_bar*np.dot((samp_size-1.0)*LD, b_bar))
	elbo_term4 = - np.sum(np.dot(bb, (samp_size-1.0)*LD)*bb)
	elbo_term5 = np.sum(np.dot(np.diag(LD*(samp_size-1)), np.transpose(postb2)))
	elbo_term7 = (-samp_size/2.0)*np.log(2.0*np.pi)

	elbo = elbo_term7 - .5*(elbo_term1 + elbo_term2 + elbo_term3 + elbo_term4 + elbo_term5) - np.sum(KL_terms)
	return elbo

def fill_in_causal_effect_size_matrix(null_mat, bs_eqtls_pmces_sparse):
	null_mat[bs_eqtls_pmces_sparse[:,0].astype(int), bs_eqtls_pmces_sparse[:,1].astype(int)] = bs_eqtls_pmces_sparse[:,2]
	return null_mat

def extract_susie_obj_from_this_bootstrap(alpha_phis, beta_phis, bs_iter, num_genes, num_snps, num_components):
	final = np.zeros((num_components, num_genes+num_snps))
	for l_iter in range(num_components):
		final[l_iter,:] = np.hstack((alpha_phis[l_iter][bs_iter,:], beta_phis[l_iter][bs_iter,:]))
	return final

def extract_susie_kl_terms_from_this_bootstrap(KL_terms, bs_iter, num_components):
	final = []
	for l_iter in range(num_components):
		final.append(KL_terms[l_iter][bs_iter])
	return np.asarray(final)


def merge_two_bootstrapped_tgfms_based_on_elbo(tgfm_obj, tgfm_obj2, variant_z_vec, variant_ld_mat, gwas_sample_size):
	# Perform seperately for each posterior sample
	# Loop through posterior samples
	for bs_iter in range(tgfm_obj.n_bs):

		# Extract causal eqtl effects for this bootstrap
		bs_eqtls_pmces = np.zeros((tgfm_obj.G, tgfm_obj.K))
		bs_eqtls_pmces_sparse = tgfm_obj.sparse_sampled_gene_eqtl_pmces[bs_iter]
		bs_eqtls_pmces = fill_in_causal_effect_size_matrix(bs_eqtls_pmces, bs_eqtls_pmces_sparse)
		# Get z-score of gene-tissue pairs and non-mediated vairants
		gene_variant_z_vec = np.hstack((tgfm_obj.nominal_twas_z[bs_iter], tgfm_obj.gwas_variant_z))
		# Extract full ld between genes, variants, and gene-variants
		gene_variant_full_ld = extract_full_gene_variant_ld(bs_eqtls_pmces, tgfm_data['reference_ld'])

		# Compute ELBO for object 1
		# Extract standard susie matrices for this posterior sample/bootrap
		bs_phi = extract_susie_obj_from_this_bootstrap(tgfm_obj.alpha_phis, tgfm_obj.beta_phis, bs_iter, tgfm_obj.G, tgfm_obj.K, tgfm_obj.L)
		bs_mu = extract_susie_obj_from_this_bootstrap(tgfm_obj.alpha_mus, tgfm_obj.beta_mus, bs_iter, tgfm_obj.G, tgfm_obj.K, tgfm_obj.L)
		bs_var = extract_susie_obj_from_this_bootstrap(tgfm_obj.alpha_vars, tgfm_obj.beta_vars, bs_iter, tgfm_obj.G, tgfm_obj.K, tgfm_obj.L)
		bs_mu2 = np.square(bs_mu) + bs_var
		bs_kl = extract_susie_kl_terms_from_this_bootstrap(tgfm_obj.KL_terms, bs_iter, tgfm_obj.L)
		# Calculate the elbo
		elbo = elbo_calc(gene_variant_z_vec, gene_variant_full_ld, gwas_sample_size, bs_phi, bs_mu, bs_mu2, bs_kl)

		# Compute ELBO for object 1
		# Extract standard susie matrices for this posterior sample/bootrap
		bs_phi_2 = extract_susie_obj_from_this_bootstrap(tgfm_obj2.alpha_phis, tgfm_obj2.beta_phis, bs_iter, tgfm_obj2.G, tgfm_obj2.K, tgfm_obj2.L)
		bs_mu_2 = extract_susie_obj_from_this_bootstrap(tgfm_obj2.alpha_mus, tgfm_obj2.beta_mus, bs_iter, tgfm_obj2.G, tgfm_obj2.K, tgfm_obj2.L)
		bs_var_2 = extract_susie_obj_from_this_bootstrap(tgfm_obj2.alpha_vars, tgfm_obj2.beta_vars, bs_iter, tgfm_obj2.G, tgfm_obj2.K, tgfm_obj2.L)
		bs_mu2_2 = np.square(bs_mu_2) + bs_var_2
		bs_kl_2 = extract_susie_kl_terms_from_this_bootstrap(tgfm_obj2.KL_terms, bs_iter, tgfm_obj2.L)
		# Calculate the elbo
		elbo_2 = elbo_calc(gene_variant_z_vec, gene_variant_full_ld, gwas_sample_size, bs_phi_2, bs_mu_2, bs_mu2_2, bs_kl_2)

		# Only need to update tgfm_obj if this is the case
		if elbo_2 > elbo:
			for l_iter in range(tgfm_obj.L):
				# Update Phi
				tgfm_obj.alpha_phis[l_iter][bs_iter, :] = np.copy(tgfm_obj2.alpha_phis[l_iter][bs_iter, :])
				tgfm_obj.beta_phis[l_iter][bs_iter, :] = np.copy(tgfm_obj2.beta_phis[l_iter][bs_iter, :])
				# Update mu
				tgfm_obj.alpha_mus[l_iter][bs_iter, :] = np.copy(tgfm_obj2.alpha_mus[l_iter][bs_iter, :])
				tgfm_obj.beta_mus[l_iter][bs_iter, :] = np.copy(tgfm_obj2.beta_mus[l_iter][bs_iter, :])
				# Update mu_var
				tgfm_obj.alpha_vars[l_iter][bs_iter, :] = np.copy(tgfm_obj2.alpha_vars[l_iter][bs_iter, :])
				tgfm_obj.beta_vars[l_iter][bs_iter, :] = np.copy(tgfm_obj2.beta_vars[l_iter][bs_iter, :])
				# Update LBF
				tgfm_obj.alpha_lbfs[l_iter][bs_iter, :] = np.copy(tgfm_obj2.alpha_lbfs[l_iter][bs_iter, :])
				tgfm_obj.beta_lbfs[l_iter][bs_iter, :] = np.copy(tgfm_obj2.beta_lbfs[l_iter][bs_iter, :])

	# Recompute PIPs
	tgfm_obj.compute_pips()

	return tgfm_obj




def tgfm_inference_shell(tgfm_data, gene_log_prior, var_log_prior, n_components, gene_tissue_pip_threshold):
	# Vector of marginal z-scores across variants
	variant_z = tgfm_data['gwas_beta']/tgfm_data['gwas_beta_se']

	# Run tgfm (with null/default initialization)
	tgfm_obj = tgfm.TGFM(L=n_components, gene_log_pi=gene_log_prior, variant_log_pi=var_log_prior, max_iter=5)
	tgfm_obj.fit(tgfm_data_obj=tgfm_data)

	# If we have a potential high PIP gene-tissue pair, we perform a second run of TGFM with an alternative initialization
	# For each posterior sample, we then use the TGFM with the larger elbo
	# The alternative initialization involves running SuSiE fine-mapping with only non-mediated variants.
	# And then using the output of that fine-mapping as initialization for the second round of TGFM. 
	# In practice, this reduced false positives in simulations (presumably caused by incomplete maximization of the variational objective)
	if np.max(np.max(tgfm_obj.expected_alpha_pips)) > gene_tissue_pip_threshold:

		# Run SuSiE with only variants to create initialization for second run of tgfm
		susie_variant_only = susieR_pkg.susie_rss(z=variant_z.reshape((len(variant_z),1)), R=tgfm_data['reference_ld'], n=tgfm_data['gwas_sample_size'], L=n_components, estimate_residual_variance=False)

		# Number of snps and number of genes in window
		num_snps = len(variant_z)
		num_genes = len(tgfm_data['genes'])

		# Create initial values of alpha, mu, and mu_var based on SuSiE-variant only analysis
		# This works because variants are constant across all bootrapped samples
		# Alpha
		alpha_init = np.zeros((n_components, num_genes + num_snps))
		alpha_init[:, num_genes:] = susie_variant_only.rx2('alpha')
		# Mu
		mu_init = np.zeros((n_components, num_genes + num_snps))
		mu_init[:, num_genes:] = susie_variant_only.rx2('mu')
		# Mu-var
		mu_var_init = np.zeros((n_components, num_genes + num_snps))
		variant_only_mu_var = susie_variant_only.rx2('mu2') - np.square(susie_variant_only.rx2('mu'))
		mu_var_init[:, num_genes:] = variant_only_mu_var

		# Delete unneccessary data from memory
		del susie_variant_only

		# Run tgfm sampler with variant init
		tgfm_obj_variant_init = tgfm.TGFM(L=n_components, gene_log_pi=gene_log_prior, variant_log_pi=var_log_prior, max_iter=5)
		tgfm_obj_variant_init.fit(tgfm_data_obj=tgfm_data, phi_init=alpha_init, mu_init=mu_init, mu_var_init=mu_var_init)

		# For each posterior sample, select the TGFM run with the larger elbo
		tgfm_obj = merge_two_bootstrapped_tgfms_based_on_elbo(tgfm_obj, tgfm_obj_variant_init, variant_z, tgfm_data['reference_ld'], tgfm_data['gwas_sample_size'])
		
	return tgfm_obj






def extract_tissue_specific_log_prior_probabilities_window_genetic_elements(tissue_specific_prior, variant_names, tissue_names):
	# Number of variants in window
	n_var = len(variant_names)
	# Number of genes in window
	n_genes = len(tissue_names)
	# Number of bootstraps done
	n_bs = len(tissue_specific_prior['variant'])

	# Initialize matrix to keep track of probabilities of each genetic element
	var_probs = np.zeros((n_var, n_bs))
	gene_probs = np.zeros((n_genes, n_bs))

	# Loop though variants
	for var_iter in range(n_var):
		var_probs[var_iter, :] = tissue_specific_prior['variant']
	# Loop through gene-tissue pairs
	for gene_iter, tissue_name in enumerate(tissue_names):
		gene_probs[gene_iter, :] = tissue_specific_prior[tissue_name]

	# Normalize rows (across both variants and gene-tissue pairs)
	normalizers = np.sum(gene_probs,axis=0) + np.sum(var_probs,axis=0)
	norm_var_probs = var_probs/normalizers
	norm_gene_probs = gene_probs/normalizers

	return np.log(norm_var_probs), np.log(norm_gene_probs)



def component_in_middle_of_window(alpha_phi_vec, beta_phi_vec, middle_gene_arr, middle_variant_arr):
	middle_gene_indices_dicti = {}
	for middle_gene_index in middle_gene_arr:
		middle_gene_indices_dicti[middle_gene_index] = 1
	middle_variant_indices_dicti = {}
	for middle_variant_index in middle_variant_arr:
		middle_variant_indices_dicti[middle_variant_index] = 1

	booler = False
	# Gene wins
	if np.max(alpha_phi_vec) > np.max(beta_phi_vec):
		best_index = np.argmax(alpha_phi_vec)
		if best_index in middle_gene_indices_dicti:
			booler = True
	else:
		best_index = np.argmax(beta_phi_vec)
		if best_index in middle_variant_indices_dicti:
			booler = True
	return booler



def load_in_tissue_specific_prior(tissue_specific_prior_summary_file):
	# Initialize dictionary to map from genetic element class name to prior probability distribution
	element_class_to_prior_probs = {}

	# Stream prior probability file
	f = open(tissue_specific_prior_summary_file)
	head_count = 0  # Used to skip header
	for line in f:

		# Parse line
		line = line.rstrip()
		data = line.split('\t')

		# Skip header
		if head_count == 0:
			head_count = head_count +1
			continue

		# class of genetic element
		class_name = data[0]
		# And corresponding bootstrapped (bs) probability distribution
		bs_probs = np.asarray(data[2].split(';')).astype(float)

		# Set probs of 0 to 1e-100
		new_bs_probs = []
		for bs_prob in bs_probs:
			if bs_prob == 0:
				new_bs_probs.append(1e-100)
			else:
				new_bs_probs.append(bs_prob)
		# Convert python array to numpy list
		new_bs_probs = np.asarray(new_bs_probs)

		# Add mapping
		element_class_to_prior_probs[class_name] = new_bs_probs
	f.close()

	return element_class_to_prior_probs


def filter_pips_to_middle_genetic_elements(genetic_element_pips, genetic_element_names, middle_genetic_elements, pip_threshold=0.01):
	middle_pips = []
	middle_names = []
	for genetic_element_iter, genetic_element_name in enumerate(genetic_element_names):
		if genetic_element_name not in middle_genetic_elements:
			continue
		if genetic_element_pips[genetic_element_iter] < pip_threshold:
			continue
		middle_pips.append(genetic_element_pips[genetic_element_iter])
		middle_names.append(genetic_element_name)
	middle_pips = np.asarray(middle_pips)
	middle_names = np.asarray(middle_names)

	indices = np.argsort(-middle_pips)
	ordered_middle_pips = middle_pips[indices]
	ordered_middle_names = middle_names[indices]


	return ordered_middle_pips, ordered_middle_names


def write_gene_tissue_pair_pips_to_output(t_gt, window_name, gt_pairs, gt_pips, middle_genetic_elements):
	# Loop through gene-tissue pairs
	for gt_index, gt_pair in enumerate(gt_pairs):
		
		# Get pip of this gene-tissue pair
		gt_pip = gt_pips[gt_index]

		# Throw out gene-tissue pairs not in middle of 1 MB region
		if gt_pair not in middle_genetic_elements:
			continue

		# Write gene-tissue pair pip to output
		t_gt.write(window_name + '\t' + gt_pair + '\t' + str(gt_pip) + '\n')

	return t_gt





###########################
# Parse command line args
###########################
parser = argparse.ArgumentParser()
parser.add_argument('--tgfm-input-data', default='None', type=str,
                    help='File summarizing which windows to run TGFM on')
parser.add_argument('--n-components', default=10, type=int,
                    help='Number of TGFM components')
parser.add_argument('--trait-name', default='None', type=str,
                    help='Name of GWAS trait')
parser.add_argument('--tissue-summary-file', default='None', type=str,
                    help='File containing information on each tissue and info on where that tissues gene models are stored')
parser.add_argument('--tissue-specific-prior', default='None', type=str,
                    help='File location of output of TGFM tissue-specific prior')
parser.add_argument('--p-value-threshold', default=1e-5, type=float,
                    help='Only run TGFM on windows with at least 1 variant with p-value less than this threshold')
parser.add_argument('--gene-tissue-pip-threshold', default=0.2, type=float,
                    help='Threshold used by TGFM to run TGFM with alternative initializations')
parser.add_argument('--parallel-job-identifier', default='None', type=str,
                    help='String corresponding to name of this parallel run')
parser.add_argument('--out', default='None', type=str,
                    help='Output stem to print tgfm results to')
args = parser.parse_args()


# Load in file containing TGFM input data
# Each row of the file corresponds to a genomic region/window to run TGFM on
# Columns correspond to tgfm input data files
tgfm_input_data = np.loadtxt(args.tgfm_input_data,dtype=str,delimiter='\t')
tgfm_input_data = tgfm_input_data[1:,:]  # Skip first line (which is header)

# Extract n_windows on this run
n_windows = tgfm_input_data.shape[0]

# Extract tissue information
# ordered_tissue_names is a numpy array of length number of tissues where each element is a string corresponding to the name of a tissue
# tissue_to_position_mapping is a dictionary mapping from tissue name to tissue position/index
ordered_tissue_names, tissue_to_position_mapping = extract_tissue_names(args.tissue_summary_file)

# Load in tissue-specific prior
tissue_specific_prior = load_in_tissue_specific_prior(args.tissue_specific_prior + '_tissue_specific_prior_summary.txt')

# Open PIP-summary file handle
pip_output_file = args.out + '_' + args.parallel_job_identifier + '_tgfm_pip_summary.txt'
t_pip = open(pip_output_file,'w')
t_pip.write('window_name\ttop_genetic_elements\tPIPs\n')

# Open gene-tissue summary file handle
gt_output_file = args.out + '_' + args.parallel_job_identifier + '_gene_tissue_pair_tgfm_pip_summary.txt'
t_gt = open(gt_output_file,'w')
t_gt.write('window_name\tgene-tissue-pair\tPIP\n')


# Loop through windows
for window_iter in range(n_windows):

	##############################
	# Extract relevent data corresponding to this window
	###############################
	data = tgfm_input_data[window_iter, :]
	# Name of window
	window_name = data[0]
	# LD file
	ld_file = data[1]
	# TGFM input data
	tgfm_input_pkl = data[2]
	# TGFM gwas input data
	tgfm_gwas_input_pkl = data[3]

	# Load in tgfm input data
	g = open(tgfm_input_pkl, "rb")
	tgfm_data = pickle.load(g)
	g.close()
	# Load in tgfm trait-gwas input data
	g = open(tgfm_gwas_input_pkl, "rb")
	tgfm_gwas_data = pickle.load(g)
	g.close()

	# get index corresponding to trait
	trait_index = np.where(tgfm_gwas_data['gwas_study_names'] == args.trait_name)[0][0]
	# Extract trait info for this trait
	gwas_sample_size = tgfm_gwas_data['gwas_sample_size'][trait_index]
	gwas_beta = tgfm_gwas_data['gwas_beta'][trait_index,:]
	gwas_beta_se = tgfm_gwas_data['gwas_beta_se'][trait_index,:]

	# Extract gwas p
	gwas_z = gwas_beta/gwas_beta_se
	gwas_p = scipy.stats.norm.sf(abs(gwas_z))*2.0

	# Ignore windows with no pvalues less than some threshold
	if np.min(gwas_p) > args.p_value_threshold:
		print('window skipped because of window pvalue threshold')
		t_pip.write(window_name + '\tNA\tNA\n')
		continue

	# Add gwas_beta to tgfm_data obj
	tgfm_data['gwas_beta'] = gwas_beta
	tgfm_data['gwas_beta_se'] = gwas_beta_se
	tgfm_data['gwas_sample_size'] = gwas_sample_size
	del tgfm_gwas_data


	# Load in variant LD
	tgfm_data['reference_ld'] = np.load(ld_file)

	# Extract log prior probabilities for genetic elements according to tissue-specific prior
	var_log_prior, gene_log_prior = extract_tissue_specific_log_prior_probabilities_window_genetic_elements(tissue_specific_prior, tgfm_data['variants'], tgfm_data['tissues'])


	##############################
	# Run TGFM
	###############################
	tgfm_obj = tgfm_inference_shell(tgfm_data, gene_log_prior, var_log_prior, args.n_components, args.gene_tissue_pip_threshold)

	##############################
	# Organize TGFM data and print to results
	###############################
	# Extract names of genetic elements
	genetic_element_names = np.hstack((tgfm_data['gene_tissue_pairs'], tgfm_data['variants']))
	# Extract dictionary list of genetic elements in the middel of this window
	middle_genetic_elements = extract_middle_genetic_elements(tgfm_data['gene_tissue_pairs'], tgfm_data['middle_gene_indices'], tgfm_data['variants'], tgfm_data['middle_variant_indices'])
	# Extract genetic element pips
	genetic_element_pips = np.hstack((tgfm_obj.expected_alpha_pips, tgfm_obj.expected_beta_pips))

	# Extract genetic elements and pips only corresponding to middle genetic elements
	# Filter at those genetic elements with pip less than pip threshold
	ordered_middle_pips, ordered_middle_genetic_elements = filter_pips_to_middle_genetic_elements(genetic_element_pips, genetic_element_names, middle_genetic_elements, pip_threshold=0.01)
	
	# Write high pip genetic elements output
	if len(ordered_middle_pips) > 0:
		# If there exist genetic elements with pip greater than or equal to threshold 0.01
		t_pip.write(window_name + '\t')
		t_pip.write(';'.join(ordered_middle_genetic_elements) + '\t')
		t_pip.write(';'.join(ordered_middle_pips.astype(str)) + '\n')
	else:
		# If there don't exist genetic elements with pip greater than or equal to threshold 0.01
		t_pip.write(window_name + '\t' + 'no_genetic_elements_pass_pip_threshold\tNA\n')
	t_pip.flush()	# Write to credible set output

	# Write gene-tissue pairs to output
	t_gt = write_gene_tissue_pair_pips_to_output(t_gt, window_name, tgfm_data['gene_tissue_pairs'], tgfm_obj.expected_alpha_pips, middle_genetic_elements)
	t_gt.flush()

	# Save all TGFM results to pkl
	tgfm_results = {}
	tgfm_results['variants'] = tgfm_data['variants']
	tgfm_results['gene_tissue_pairs'] = tgfm_data['gene_tissue_pairs']
	tgfm_results['genes'] = tgfm_data['genes']
	tgfm_results['tissues'] = tgfm_data['tissues']
	tgfm_results['alpha_phis'] = tgfm_obj.alpha_phis
	tgfm_results['alpha_mus'] = tgfm_obj.alpha_mus
	tgfm_results['alpha_vars'] = tgfm_obj.alpha_vars
	tgfm_results['beta_phis'] = tgfm_obj.beta_phis
	tgfm_results['expected_alpha_pips'] = tgfm_obj.expected_alpha_pips
	tgfm_results['expected_beta_pips'] = tgfm_obj.expected_beta_pips
	tgfm_results['nominal_twas_z'] = tgfm_obj.nominal_twas_z
	tgfm_results['middle_variant_indices'] = tgfm_data['middle_variant_indices']
	tgfm_results['middle_gene_indices'] = tgfm_data['middle_gene_indices']
	

	# Write pickle file
	window_tgfm_output_file = args.out + '_' + window_name + '_tgfm_results.pkl'
	g = open(window_tgfm_output_file, "wb")
	pickle.dump(tgfm_results, g)
	g.close()



# Close file handles
t_pip.close()
t_gt.close()