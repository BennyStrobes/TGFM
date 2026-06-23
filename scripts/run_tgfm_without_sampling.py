import sys
import pandas as pd
import numpy as np 
import os 
import pdb
import pickle
import scipy.stats
import rpy2
import rpy2.robjects.numpy2ri as numpy2ri
import rpy2.robjects as ro
ro.conversion.py2ri = numpy2ri
numpy2ri.activate()
from rpy2.robjects.packages import importr
susieR_pkg = importr('susieR')
import argparse
import tgfm_init

def get_credible_set_genes(phi_comp, cs_thresh):
	ordered_genes = np.argsort(-phi_comp)
	cs_genes = []
	cs_counter = 0.0
	for gene_index in ordered_genes:
		if cs_counter < cs_thresh:
			cs_genes.append(gene_index)
		cs_counter = cs_counter + phi_comp[gene_index]
	cs_genes = np.asarray(cs_genes)
	return cs_genes


def extract_valid_joint_susie_components_from_full_ld(alpha_phi, beta_phi, full_ld, ld_thresh, subset_n = 100):
	num_components = alpha_phi.shape[0]
	valid_components = []

	num_genes = alpha_phi.shape[1]
	num_variants = beta_phi.shape[1]

	for component_num in range(num_components):
		cs_predictors = get_credible_set_genes(np.hstack((alpha_phi[component_num,:], beta_phi[component_num,:])), .95)

		if subset_n > len(cs_predictors):
			# absolute ld among genes and variants in credible set
			if np.min(np.abs(full_ld[cs_predictors,:][:, cs_predictors])) > ld_thresh:
				valid_components.append(component_num)
		else:
			# First run subsetted analysis
			cs_predictors_subset = np.random.choice(cs_predictors, size=subset_n, replace=False, p=None)
			if np.min(np.abs(full_ld[cs_predictors_subset,:][:, cs_predictors_subset])) > ld_thresh:
				if np.min(np.abs(full_ld[cs_predictors,:][:, cs_predictors])) > ld_thresh:
					valid_components.append(component_num)

	return valid_components





def extract_uniform_log_prior_probabilities(variant_names, gene_names):
	n_snps = len(variant_names)
	n_genes = len(gene_names)
	n_ele = n_snps + n_genes

	ele_prob = 1.0/n_ele

	snp_prior = np.ones(n_snps)*ele_prob
	gene_prior = np.ones(n_genes)*ele_prob

	return np.log(snp_prior), np.log(gene_prior)


def extract_full_gene_variant_ld(standardized_eqtl_effects, variant_ld):
	expression_covariance = np.dot(np.dot(standardized_eqtl_effects, variant_ld), np.transpose(standardized_eqtl_effects))
	np.fill_diagonal(expression_covariance, 1.0)
	dd = np.diag(1.0/np.sqrt(np.diag(expression_covariance)))
	ge_ld = np.dot(np.dot(dd, expression_covariance),dd)
	gene_variant_ld = np.dot(standardized_eqtl_effects,variant_ld) # Ngenes X n_variants
	top = np.hstack((ge_ld, gene_variant_ld))
	bottom = np.hstack((np.transpose(gene_variant_ld), variant_ld))
	full_ld = np.vstack((top,bottom))
	return full_ld

def initialize_tgfm_object(tgfm_data, gene_log_prior, var_log_prior, n_components):
	# Initialize tgfm results object
	# Currently a bit hacky
	# All this does is make a python data object (useful for later)
	tgfm_data_obj_light = {}
	tgfm_data_obj_light['gene_tissue_pairs'] = np.copy(tgfm_data['gene_tissue_pairs'])
	tgfm_data_obj_light['variants'] = np.copy(tgfm_data['variants'])
	tgfm_data_obj_light['gwas_sample_size'] = np.copy(tgfm_data['gwas_sample_size'])
	tgfm_obj = tgfm_init.TGFM(L=n_components, gene_init_log_pi=gene_log_prior, variant_init_log_pi=var_log_prior)
	tgfm_obj.fit(twas_data_obj=tgfm_data_obj_light)

	return tgfm_obj

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

	# Compute pips
	pips = compute_pips(susie_obj.rx2('alpha'))
	tgfm_obj.alpha_pip = pips[:(tgfm_obj.G)]
	tgfm_obj.beta_pip = pips[(tgfm_obj.G):]

	# Log bayes factors
	lbf = susie_obj.rx2('lbf_variable')
	tgfm_obj.alpha_lbf = lbf[:, :(tgfm_obj.G)]
	tgfm_obj.beta_lbf = lbf[:, (tgfm_obj.G):]

	return tgfm_obj


def compute_pips(alpha_mat):
	LL = alpha_mat.shape[0]
	n_elements = alpha_mat.shape[1]
	anti_pips = np.ones(n_elements)

	for component_iter in range(LL):
		anti_pips = anti_pips*(1.0 - alpha_mat[component_iter,:])
	pips = 1.0 - anti_pips
	return pips


def tgfm_no_sampling_inference_shell(tgfm_data, gene_log_prior, var_log_prior, gene_variant_full_ld, n_components, gene_tissue_pip_threshold):
	# Initialize TGFM object
	tgfm_obj = initialize_tgfm_object(tgfm_data, gene_log_prior, var_log_prior, n_components)

	# Get z-scores of each genetic element
	variant_z = tgfm_data['gwas_beta']/tgfm_data['gwas_beta_se']
	new_gene_z = np.dot(tgfm_data['gene_eqtl_pmces'], variant_z)  # No need to normalize by genetic variance of gene because eqtl_pmces are standardized to have variance 1
	# Add gene z-scores to tgfm data object
	tgfm_obj.nominal_twas_z = new_gene_z
	# Create vector of concatenated z-scores
	z_vec = np.hstack((new_gene_z,variant_z))

	# Create concatenated vector of prior probs
	prior_probs = np.hstack((np.exp(gene_log_prior), np.exp(var_log_prior)))

	# Add warnings
	rpy2.robjects.r['options'](warn=1)

	# Run susie jointly across gene-tissue pairs and non-mediated variants (with null/default initialization)
	susie_null_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=n_components, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=False)
	susie_null_init_elbo = susie_null_init.rx2('elbo')[-1]
	# Add susie fine-mapping results to tgfm data object
	tgfm_obj = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_null_init)
	# Delete unnecessary data from memory
	del susie_null_init

	# If we have a potential high PIP gene-tissue pair, we perform a second run of TGFM with an alternative initialization
	# We then use the TGFM with the larger elbo
	# The alternative initialization involves running SuSiE fine-mapping with only non-mediated variants.
	# And then using the output of that fine-mapping as initialization for the second round of TGFM. 
	# In practice, this reduced false positives in simulations (presumably caused by incomplete maximization of the variational objective)
	if np.max(tgfm_obj.alpha_pip) > gene_tissue_pip_threshold:
		
		# Compute prior for variant only analysis
		# basically equal weight on variants and 0 weight on gene-tissue pairs
		p_var_only = np.ones(len(z_vec))
		p_var_only[:len(tgfm_obj.nominal_twas_z)] = 0.0
		p_var_only = p_var_only/np.sum(p_var_only)
		# Run susie with only variants
		susie_variant_only = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=n_components, prior_weights=p_var_only.reshape((len(p_var_only),1)), estimate_residual_variance=False)

		# Second run of TGFM with variant initialization
		init_obj = {'alpha':susie_variant_only.rx2('alpha'), 'mu':susie_variant_only.rx2('mu'),'mu2':susie_variant_only.rx2('mu2')}
		init_obj2 = ro.ListVector(init_obj)
		# Delete unneccessary data from memory
		del susie_variant_only
		init_obj2.rclass = rpy2.robjects.StrVector(("list", "susie"))
		susie_variant_init = susieR_pkg.susie_rss(z=z_vec.reshape((len(z_vec),1)), R=gene_variant_full_ld, n=tgfm_data['gwas_sample_size'], L=n_components, s_init=init_obj2, prior_weights=prior_probs.reshape((len(prior_probs),1)), estimate_residual_variance=False)

		# Get elbo of the second run of tgfm
		susie_variant_init_elbo = susie_variant_init.rx2('elbo')[-1]

		# Select model with largest elbo
		if susie_variant_init_elbo > susie_null_init_elbo:
			# Variant init wins (don't need to update if the null init wins)
			tgfm_obj = update_tgfm_obj_with_susie_res_obj(tgfm_obj, susie_variant_init)
		# Delete unneccessary data from memory
		del susie_variant_init
	
	return tgfm_obj

def extract_middle_genetic_elements(ordered_genes, middle_gene_indices, ordered_variants, middle_variant_indices):
	dicti = {}
	for gene_name in ordered_genes[middle_gene_indices]:
		dicti[gene_name] = 1
	for variant_name in ordered_variants[middle_variant_indices.astype(int)]:
		dicti[variant_name] = 1
	return dicti

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
tgfm_input_data = tgfm_input_data[1:,:]

# Get n_windows on this run
n_windows = tgfm_input_data.shape[0]



# Open PIP file handle
pip_output_file = args.out + '_' + args.parallel_job_identifier + '_tgfm_no_sample_pip_summary.txt'
t_pip = open(pip_output_file,'w')
t_pip.write('window_name\ttop_genetic_elements\tPIPs\n')



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

	# Load in LD
	ld_mat = np.load(ld_file)

	# Extract full ld between genes, variants, and gene-variants
	gene_variant_full_ld = extract_full_gene_variant_ld(tgfm_data['gene_eqtl_pmces'], ld_mat)
	del ld_mat

	# Extract prior probabilities of each genetic element being fine-mapped
	# Here, we use a uniform prior
	var_log_prior, gene_log_prior = extract_uniform_log_prior_probabilities(tgfm_data['variants'], tgfm_data['gene_tissue_pairs'])


	##############################
	# Run TGFM
	###############################
	tgfm_obj = tgfm_no_sampling_inference_shell(tgfm_data, gene_log_prior, var_log_prior, gene_variant_full_ld, args.n_components, args.gene_tissue_pip_threshold)


	##############################
	# Organize TGFM data and print to results
	###############################
	# Extract TGFM components that pass purity filter
	valid_tgfm_components = extract_valid_joint_susie_components_from_full_ld(tgfm_obj.alpha_phi, tgfm_obj.beta_phi, gene_variant_full_ld, .5)

	# Extract names of genetic elements
	genetic_element_names = np.hstack((tgfm_data['gene_tissue_pairs'], tgfm_data['variants']))
	# Extract dictionary list of genetic elements in the middel of this window
	middle_genetic_elements = extract_middle_genetic_elements(tgfm_data['gene_tissue_pairs'], tgfm_data['middle_gene_indices'], tgfm_data['variants'], tgfm_data['middle_variant_indices'])


	# Extract genetic element pips
	genetic_element_pips = np.hstack((tgfm_obj.alpha_pip, tgfm_obj.beta_pip))

	# Extract genetic elements and pips only corresponding to middle genetic elements
	# Filter at those genetic elements with pip less than pip threshold
	ordered_middle_pips, ordered_middle_genetic_elements = filter_pips_to_middle_genetic_elements(genetic_element_pips, genetic_element_names, middle_genetic_elements, pip_threshold=0.01)
	
	# Write to credible set output
	if len(ordered_middle_pips) > 0:
		# If there exist genetic elements with pip greater than or equal to threshold 0.01
		t_pip.write(window_name + '\t')
		t_pip.write(';'.join(ordered_middle_genetic_elements) + '\t')
		t_pip.write(';'.join(ordered_middle_pips.astype(str)) + '\n')
	else:
		# If there don't exist genetic elements with pip greater than or equal to threshold 0.01
		t_pip.write(window_name + '\t' + 'no_genetic_elements_pass_pip_threshold\tNA\n')
	t_pip.flush()


	# Save all TGFM results to pkl
	tgfm_results = {}
	tgfm_results['variants'] = tgfm_data['variants']
	tgfm_results['gene_tissue_pairs'] = tgfm_data['gene_tissue_pairs']
	tgfm_results['genes'] = tgfm_data['genes']
	tgfm_results['tissues'] = tgfm_data['tissues']
	tgfm_results['alpha_phi'] = tgfm_obj.alpha_phi
	tgfm_results['beta_phi'] = tgfm_obj.beta_phi
	tgfm_results['alpha_lbf'] = tgfm_obj.alpha_lbf
	tgfm_results['beta_lbf'] = tgfm_obj.beta_lbf
	tgfm_results['alpha_mu'] = tgfm_obj.alpha_mu
	tgfm_results['beta_mu'] = tgfm_obj.beta_mu
	tgfm_results['alpha_var'] = tgfm_obj.alpha_var
	tgfm_results['beta_var'] = tgfm_obj.beta_var
	tgfm_results['alpha_pip'] = tgfm_obj.alpha_pip
	tgfm_results['beta_pip'] = tgfm_obj.beta_pip
	tgfm_results['component_variances'] = tgfm_obj.component_variances
	tgfm_results['valid_components'] = valid_tgfm_components
	tgfm_results['nominal_twas_z'] = tgfm_obj.nominal_twas_z
	tgfm_results['middle_variant_indices'] = tgfm_data['middle_variant_indices']
	tgfm_results['middle_gene_indices'] = tgfm_data['middle_gene_indices']

	# Write pickle file
	window_tgfm_output_file = args.out + '_' + window_name + '_tgfm_no_sample_results.pkl'
	g = open(window_tgfm_output_file, "wb")
	pickle.dump(tgfm_results, g)
	g.close()


# Close file handles
t_pip.close()
