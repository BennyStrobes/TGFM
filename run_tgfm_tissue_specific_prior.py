import sys
sys.path.remove('/n/app/python/3.7.4-ext/lib/python3.7/site-packages')
import pandas as pd
import numpy as np 
import os 
import pdb
import scipy.special
import pickle
import argparse



##**
def concatenate_tgfm_no_sampling_pip_summary_files_across_parallel_jobs(tgfm_input_data_identifiers, tgfm_without_sampling_output, concat_output_file):
	# Open Output file handle
	t = open(concat_output_file,'w')

	# Loop through results files
	for ii, tgfm_input_data_identifier in enumerate(tgfm_input_data_identifiers):
		# Results file corresponding to a single run
		job_input_file = tgfm_without_sampling_output + '_' + tgfm_input_data_identifier + '_tgfm_no_sample_pip_summary.txt'
		f = open(job_input_file)
		head_count = 0
		# Add all lines of results file to concatenated
		for line in f:
			line = line.rstrip()
			# if first parallelized result file, print header
			if head_count == 0:
				head_count = head_count + 1
				if ii == 0:
					t.write(line + '\n')
				continue
			t.write(line + '\n')
		f.close()
	t.close()
	return

## **
def component_in_middle_of_window(alpha_phi_vec, beta_phi_vec, middle_gene_indices_dicti, middle_variant_indices_dicti):
	booler = False
	# To check wither component is in the middle of the window, uses the highest propbablity genetic element as the marker
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

##***
def extract_tissue_names(tissue_summary_file):
	f = open(tissue_summary_file)
	head_count = 0
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		tissue_name = data[0]
		arr.append(tissue_name)
	f.close()
	return np.asarray(arr)

##**
def generate_component_level_summary_data(concatenated_tgfm_no_sampling_pip_summary_file, component_level_summary_file, tissue_names,tgfm_without_sampling_output, per_window_abf_output_stem, max_components_per_window):
	
	# Create mapping from tissue name to position/index
	tiss_to_position_mapping = {}
	for ii,tissue_name in enumerate(tissue_names):
		tiss_to_position_mapping[tissue_name] = ii

	# Open output file handle
	t = open(component_level_summary_file,'w')
	t.write('window_name\tcomponent_number\tgenetic_element_names\tgenetic_element_class_names\tmiddle_genetic_element_names\twindow_abf_file\n')

	# Open input file handle containing results of TGFM no sampling
	# One line for each genomic window
	f = open(concatenated_tgfm_no_sampling_pip_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')

		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue

		# Name of window
		window_name = data[0]

		# Skip windows that did not pass filters (e.g., no snps with marginal gwas pvalue less than some threshold.)
		if data[1] == 'NA':
			continue

		# Load in TGFM results pkl file for this window
		window_pkl_file = tgfm_without_sampling_output + '_' + window_name + '_tgfm_no_sample_results.pkl'
		g = open(window_pkl_file, "rb")
		tgfm_results = pickle.load(g)
		g.close()

		# Create dictionary list of indices of:
		# (1) Genes in the middle MB of the window
		# (2) Variants in the middle MB of the window
		# As only these genetic elements count towards prior
		middle_gene_indices_dicti = {}
		middle_variant_indices_dicti = {}
		for indexer in tgfm_results['middle_gene_indices']:
			middle_gene_indices_dicti[indexer] =1
		for indexer in tgfm_results['middle_variant_indices']:
			middle_variant_indices_dicti[indexer] = 1

		# Filter out windows with no identified TGFM components
		if 'valid_components' not in tgfm_results:
			continue
		# If specified, filter out windows with greater than args.max-components-per-window
		if max_components_per_window is not None and len(tgfm_results['valid_components']) > max_components_per_window:
			continue

		# Initialize potential output data for this window
		# File containing abf data for this window
		per_window_abf_output_file = per_window_abf_output_stem + window_name + '.npy'
		# And array to keep track of abf data
		per_window_abf = []
		# Boolean on wheter the window has at least 1 component in the middle MB
		window_has_component = False
		# Array to keep track of components
		window_component_lines = []

		# Loop through all tgfm components (though not all of these may pass purity filter)
		for component_iter in range(tgfm_results['alpha_lbf'].shape[0]):
			# Extract log approximate bayes factor (factor)
			log_abf = np.hstack((tgfm_results['alpha_lbf'][component_iter, :], tgfm_results['beta_lbf'][component_iter, :]))
			# Extract names of genetic elements in the window
			element_names = np.hstack((tgfm_results['gene_tissue_pairs'], tgfm_results['variants']))
			# Extract the classes of each genetic element
			class_names = np.hstack((tgfm_results['tissues'], np.asarray(['variant']*len(tgfm_results['variants']))))
			# Extract names of genetic elements in the middle of the window
			middle_element_names = np.hstack((tgfm_results['gene_tissue_pairs'][tgfm_results['middle_gene_indices']], tgfm_results['variants'][tgfm_results['middle_variant_indices']]))
			# And line to window array
			liner = window_name + '\t' + str(component_iter) + '\t' + ';'.join(element_names) + '\t' + ';'.join(class_names) + '\t' + ';'.join(middle_element_names) + '\t' + per_window_abf_output_file + '\n'
			window_component_lines.append(liner)
			per_window_abf.append(log_abf)
			# Check if component passes purity filter and is in the middle of the window
			if component_iter in tgfm_results['valid_components'] and component_in_middle_of_window(tgfm_results['alpha_phi'][component_iter, :], tgfm_results['beta_phi'][component_iter, :], middle_gene_indices_dicti, middle_variant_indices_dicti):
				window_has_component = True

		# Convert to np array and save IF the window has at least 1 component passing the purity filter
		if window_has_component:
			for liner in window_component_lines:
				t.write(liner)
			per_window_abf = np.asarray(per_window_abf)
			np.save(per_window_abf_output_file, per_window_abf)
	t.close()
	f.close()
	return




##**
def extract_alpha_vec_given_lbf_and_prior_matrix(lbf, prior_probs):
	maxlbf = np.max(lbf)
	ww = np.exp(lbf - maxlbf)
	w_weighted = (np.transpose(prior_probs)*ww)
	alpha = np.transpose(w_weighted)/np.sum(w_weighted,axis=1)
	return alpha



##**
def create_mapping_from_window_to_class_to_indices(component_level_abf_summary_file, tissue_name_to_position, tissue_names):
	# Initialize dictionary to keep track of per window mapping from genetic element class (e.g. non-mediated variant, or Adipose_Tissue) to indices of genetic elements in that class
	dicti = {}
	dicti_middle = {}
	dicti_variant_middle = {}

	# Stream comonent level summary files
	f = open(component_level_abf_summary_file)
	head_count = 0  # Used to skip header
	for line in f:
		line = line.rstrip()
		data = line.split('\t')

		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue

		# Name of window
		window_name = data[0]

		# There can be multiple components from a single window
		# No need to learn the mapping multiple times as the mapping is the same across all components in the window
		if window_name in dicti:
			continue

		# initialize per-window sub-dictionary
		sub_dicti = {}
		sub_dicti_middle = {}

		# Add keys corresponding to all classes of genetic elements
		sub_dicti['variant'] = []
		sub_dicti_middle['variant'] = []
		for tissue_name in tissue_names:
			sub_dicti[tissue_name] = []
			sub_dicti_middle[tissue_name] = []

		# Extract element names in window and middle of window
		ele_names = np.asarray(data[2].split(';'))
		class_names = np.asarray(data[3].split(';'))
		middle_ele_names = np.asarray(data[4].split(';'))
		# Create dictionary of middle genetic elements
		middle_ele_dicti = {}
		for middle_ele_name in middle_ele_names:
			middle_ele_dicti[middle_ele_name] = 1


		# loop through elements and add to sub-dictionary
		for ele_iter, ele_name in enumerate(ele_names):
			# Get class of this genetic element
			ele_class = class_names[ele_iter]

			# Add to dictionary
			sub_dicti[ele_class].append(ele_iter)
			# If element in middle MB, add to middle dictionary
			if ele_name in middle_ele_dicti:
				sub_dicti_middle[ele_class].append(ele_iter)
		
		# Convert from python lists to np arrays
		sub_dicti['variant'] = np.asarray(sub_dicti['variant'])
		sub_dicti_middle['variant'] = np.asarray(sub_dicti_middle['variant'])
		for tissue_name in tissue_names:
			sub_dicti[tissue_name] = np.asarray(sub_dicti[tissue_name])
			sub_dicti_middle[tissue_name] = np.asarray(sub_dicti_middle[tissue_name])
		
		# Add n_elements key to sub_dicti
		sub_dicti['n_elements'] = len(ele_names)

		# add window sub-dictionary to global dictionary
		dicti[window_name] = sub_dicti
		dicti_middle[window_name] = sub_dicti_middle

	f.close()
	return dicti, dicti_middle

def randomly_sample_bootrap_window_indices(n_bootstraps, n_windows):
	bs_indices = []
	for bs_iter in range(n_bootstraps):
		indices = np.random.choice(np.arange(n_windows), size=n_windows, replace=True)
		bs_indices.append(indices)

	# Create mapping from bootstrap iteration to window indices
	bs_mapping = []
	for bs_iter in range(n_bootstraps):
		indices = bs_indices[bs_iter]
		temper = {}
		for index in indices:
			if index not in temper:
				temper[index] = 1.0
			else:
				temper[index] = temper[index] + 1.0
		bs_mapping.append(temper)
	return bs_mapping

def create_mapping_from_window_name_to_bootstrap_iterations(window_names, n_bootstraps, bootstrap_iter_to_window_indices):
	# Initialize dictionary
	window_to_bootstraps = {}

	# Loop through windows
	for window_iter, window_name in enumerate(window_names):
		# Initialize array to keep track of how many times window is in each bootstrap
		vec = np.zeros(n_bootstraps)

		# Loop through bootstrap iters
		for bs_iter in range(n_bootstraps):

			# Find out if window in bootrap
			if window_iter in bootstrap_iter_to_window_indices[bs_iter]:
				# And if so, how many times is the window in the bootsrap
				vec[bs_iter] = bootstrap_iter_to_window_indices[bs_iter][window_iter]

		window_to_bootstraps[window_name] = vec

	return window_to_bootstraps



##**
def learn_iterative_variant_gene_tissue_prior_pip_level_bootstrapped(component_level_abf_summary_file, tissue_names, per_window_abf_output_stem, max_iter=400, n_bootstraps=100):
	# Create mapping from tissue name to tissue position/index
	tissue_name_to_position = {}
	for ii, tissue_name in enumerate(tissue_names):
		tissue_name_to_position[tissue_name] = ii

	# Initialize prior probability emperical distribution
	variant_prob_distr = np.ones(n_bootstraps)*.1
	tissue_probs_distr = np.ones((len(tissue_names), n_bootstraps))*.1

	#Create mappings to keep track of per window mapping from genetic element class (e.g. non-mediated variant, or Adipose_Tissue) 
	# to indices of genetic elements in that class
	window_to_class_to_indices, window_to_class_to_middle_indices = create_mapping_from_window_to_class_to_indices(component_level_abf_summary_file, tissue_name_to_position, tissue_names)

	# Get ordered window names
	window_names = np.asarray([*window_to_class_to_indices])
	# And the number of windows
	n_windows = len(window_names)


	# Randomly sample (with replacement) which windows correspond to each bootstrap
	# Return python list of length number of bootstraps
	# Where each element of list is a dictionary list of window indexes included in that bootstrapping iterations.
	# The keys are window indices and values are the number of times that window is included in the bootstrap
	bootstrap_iter_to_window_indices = randomly_sample_bootrap_window_indices(n_bootstraps, n_windows)

	# Go in reverse direction: Create mapping from window_name to bootstrap iterations
	# Create mapping from window name to bootstrap iterations
	# Return a dictionary where keys are names of windows and values are vector of length number of bootstraps where elements tell the number of times the window is each bootrap
	window_to_bootstraps = create_mapping_from_window_name_to_bootstrap_iterations(window_names, n_bootstraps, bootstrap_iter_to_window_indices)


	# Iterative algorithm
	for itera in range(max_iter):
		variant_prob_distr, tissue_probs_distr = update_prior_prob_for_variant_gene_tissue_bootstrapped(component_level_abf_summary_file, tissue_name_to_position, variant_prob_distr, tissue_probs_distr, tissue_names, window_to_class_to_indices, window_to_class_to_middle_indices, n_bootstraps, window_names, window_to_bootstraps)


	return variant_prob_distr, tissue_probs_distr




def compute_expected_pips(alpha_mats):
	if len(alpha_mats) == 0:
		print('assumption errorr')
		pdb.set_trace()
	# Initialize expected pips
	pips = np.ones(alpha_mats[0].shape)
	for component_iter in range (len(alpha_mats)):
		pips = pips*(1.0-alpha_mats[component_iter])
	pips = 1.0 - pips

	return pips


##**
def update_prior_prob_for_variant_gene_tissue_bootstrapped(component_level_abf_summary_file, tissue_name_to_position, variant_prob_distr, tissue_probs_distr, tissue_names, window_to_class_to_indices, window_to_class_to_middle_indices, n_bootstraps, window_names, window_to_bootstraps):
	# General info
	n_tissues = len(tissue_names)	
	n_windows = len(window_names)

	# Keep track of number of genetic elements in each class for each bootstrap
	variant_counts = np.zeros(n_bootstraps)
	tissue_counts = np.zeros((n_tissues, n_bootstraps))
	# Keep track of sum of genetic element posterior probabilities in each class for each bootstrap
	variant_posterior_sum = np.zeros(n_bootstraps)
	tissue_posterior_sum = np.zeros((n_tissues, n_bootstraps))

	# Create dictionary mapping (hash) from genetic element class to prior probability of said class
	# Prior probabilities updated each iterations
	prior_hash = {}
	prior_hash['variant'] = variant_prob_distr
	for ii in range(n_tissues):
		prior_hash[tissue_names[ii]] = tissue_probs_distr[ii,:]

	# Create mapping from window name to prior probability of each genetic element
	window_to_prior_probs = {}
	# Do this for each window
	for window_name in window_names:
		# Initialize prior probability matrix
		prior_prob_raw = np.zeros((window_to_class_to_indices[window_name]['n_elements'], n_bootstraps))

		# First do variants
		# all variants have the same prior
		class_name = 'variant'
		class_indices = window_to_class_to_indices[window_name][class_name]
		prior_prob_raw[class_indices, :] = prior_hash[class_name]

		# Now add priors for each tissue
		for class_name in tissue_names:
			class_indices = window_to_class_to_indices[window_name][class_name]
			if len(class_indices) == 0:
				continue
			prior_prob_raw[class_indices, :] = prior_hash[class_name]

		# Normalize prior probabalities to sume to one within each bootstrapped sample
		prior_probs = prior_prob_raw/np.sum(prior_prob_raw,axis=0)

		# Add window prior_probs to global mapping
		window_to_prior_probs[window_name] = prior_probs


	# Note: the following code relies on the fact that all components from a single window are all next to each other in the component_level_abf_summary_file
	# Loop through components
	head_count = 0
	prev_window_name = 'NULL'
	f = open(component_level_abf_summary_file)
	for line in f:

		# parse line
		line = line.rstrip()
		data = line.split('\t')

		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue

		# Name of window
		window_name = data[0]


		# We are at a component in a NEW window: Current component is from a different window than previous component
		if window_name != prev_window_name:

			# If not the first window
			if prev_window_name != 'NULL':
			
				# Compute pips of all genetic elements in the window
				expected_pips = compute_expected_pips(alpha_mats)
	
				# Update tally of (based on new pips from this window):
				### 1. number of genetic elements in each class for each bootstrap
				### 2. sum of genetic element posterior probabilities in each class for each bootstrap
				# Variants first
				indices = np.asarray(window_to_class_to_middle_indices[prev_window_name]['variant'])  # indexes of genetic elements that are variants
				bs_scaling_factors = window_to_bootstraps[prev_window_name]  # Number of times window was in each bootstrap
				variant_posterior_sum = variant_posterior_sum + (bs_scaling_factors*np.sum(expected_pips[indices,:],axis=0))
				variant_counts = variant_counts + (bs_scaling_factors*len(indices))
				# Now genes
				for tiss_iter, tissue_name in enumerate(tissue_names):
					indices = np.asarray(window_to_class_to_middle_indices[prev_window_name][tissue_name])  # indexes of genetic elements that are genes from current tissue
					if len(indices) == 0:  # Skip tissue if no genes from that tissue
						continue
					tissue_posterior_sum[tiss_iter, :] = tissue_posterior_sum[tiss_iter, :] + (bs_scaling_factors*np.sum(expected_pips[indices,:],axis=0))
					tissue_counts[tiss_iter, :] = tissue_counts[tiss_iter, :] + (bs_scaling_factors*len(indices))

			# Load in data for new window
			window_abfs = np.load(data[5])  # Matrix of approximate bayes factors (abfs). Rows are components, columns are genetic elements
			component_counter = 0
			prev_window_name = window_name
			alpha_mats = []

		# Approximate bayes factor for current component
		ele_labfs = window_abfs[component_counter, :]

		# Keep track of which component we are at
		component_counter = component_counter + 1

		# Get prior probabilities of all genetic elements in window
		prior_probs = window_to_prior_probs[window_name]
		
		# Compute inclusion probabilities (alpha) for this component across all bootstraps
		# Bootstraps will have different values of alpha for the same component because prios are different
		# alpha_mat is matrix of dimension n_genetic_elements by num_bootstraps
		alpha_mat = extract_alpha_vec_given_lbf_and_prior_matrix(ele_labfs, prior_probs)

		# Add to arr to keep track of all alpha_mats in window (one for each component)
		alpha_mats.append(alpha_mat)
	f.close()

	# This is the same code as above following: prev_window_name != 'NULL':
	# The reason here is we need to tally results from the last window
	# Need to fix so not copied code

	# Compute pips of all genetic elements in the window
	expected_pips = compute_expected_pips(alpha_mats)

	# Update tally of (based on new pips from this window):
	### 1. number of genetic elements in each class for each bootstrap
	### 2. sum of genetic element posterior probabilities in each class for each bootstrap
	# Variants first
	indices = np.asarray(window_to_class_to_middle_indices[window_name]['variant'])  # indexes of genetic elements that are variants
	bs_scaling_factors = window_to_bootstraps[window_name]  # Number of times window was in each bootstrap
	variant_posterior_sum = variant_posterior_sum + (bs_scaling_factors*np.sum(expected_pips[indices,:],axis=0))
	variant_counts = variant_counts + (bs_scaling_factors*len(indices))
	# Now genes
	for tiss_iter, tissue_name in enumerate(tissue_names):
		indices = np.asarray(window_to_class_to_middle_indices[window_name][tissue_name])  # indexes of genetic elements that are genes from current tissue
		if len(indices) == 0:
			continue
		tissue_posterior_sum[tiss_iter, :] = tissue_posterior_sum[tiss_iter, :] + (bs_scaling_factors*np.sum(expected_pips[indices,:],axis=0))
		tissue_counts[tiss_iter, :] = tissue_counts[tiss_iter, :] + (bs_scaling_factors*len(indices))


	# Deal with fact some tissue_counts may be zero
	tissue_counts[tissue_counts == 0.0] = .1  # Avoids divide by zero problem 

	# For each bootrapped index, compute probability a variant is causal
	variant_prob_distr = variant_posterior_sum/variant_counts   # Vector of length n_bootstraps
	# Compute probability a gene from each tissue is causal
	tissue_probs_distr = tissue_posterior_sum/tissue_counts  # Matrix of dimension number of tissues by number of bootstraps

	return variant_prob_distr, tissue_probs_distr




##**
def load_in_tgfm_parallel_job_identifiers(tgfm_input_data_summary_file):
	tgfm_input_data_identifiers = []
	f = open(tgfm_input_data_summary_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		identifier = data[0]
		tgfm_input_data_identifiers.append(identifier)

	f.close()
	return np.asarray(tgfm_input_data_identifiers)




###########################
# Parse command line args
###########################
parser = argparse.ArgumentParser()
parser.add_argument('--trait-name', default='None', type=str,
                    help='Name of GWAS trait')
parser.add_argument('--tissue-summary-file', default='None', type=str,
                    help='File containing information on each tissue and info on where that tissues gene models are stored')
parser.add_argument('--tgfm-parallel-job-identifier_file', default='None', type=str,
                    help='File containing information on the names of the parallel runs of the tgfm_without_sampling analysis')
parser.add_argument('--tgfm-without-sampling-output', default='None', type=str,
                    help='Output stem of results of TGFM without sampling analysis')
parser.add_argument('--num-iterations', default=400, type=int,
                    help='Number of posterior samples to randomly draw from each gene-tissue pair')
parser.add_argument('--num-posterior-samples', default=100, type=int,
                    help='Number of posterior samples to randomly draw from each gene-tissue pair')
parser.add_argument('--max-components-per-window', default=None, type=int,
                    help='Filter out windows with more than this number of components per window. Really large numbers of components per window could be an error.')
parser.add_argument('--random-seed', default=None, type=int,
                    help='Random seed')
parser.add_argument('--out', default='None', type=str,
                    help='Output stem to print tgfm input data to')
args = parser.parse_args()



# Set random seed if requested
if args.random_seed is not None:
	np.random.seed(args.random_seed)

print('Wrangling data from TGFM without sampling analysis')

#Extract tissue names
tissue_names = extract_tissue_names(args.tissue_summary_file)


# Load in TGFM input data
tgfm_parallel_job_identifiers = load_in_tgfm_parallel_job_identifiers(args.tgfm_parallel_job_identifier_file)



# Concatenate TGFM no sampling PIP summary files across parallel runs (one line for each window)
concatenated_tgfm_no_sampling_pip_summary_file = args.out + '_tissue_specific_prior_tgfm_no_sampling_pip_summary.txt'
concatenate_tgfm_no_sampling_pip_summary_files_across_parallel_jobs(tgfm_parallel_job_identifiers, args.tgfm_without_sampling_output, concatenated_tgfm_no_sampling_pip_summary_file)


# Extract TGFM components to learn tissue-specific prior off of 
# Create a file with a line for each component
# For each component, extract the approximate bayes factors (abf) of each genetic element from the TGFM no sampling results
component_level_summary_file = args.out + '_tissue_specific_prior_abf_summary.txt'
per_window_abf_output_stem = args.out + '_tissue_specific_prior_per_window_abf_'
generate_component_level_summary_data(concatenated_tgfm_no_sampling_pip_summary_file, component_level_summary_file, tissue_names, args.tgfm_without_sampling_output,  per_window_abf_output_stem, args.max_components_per_window)



###################################################
# Learn TGFM Tissue-specific prior
###################################################
print('Running (iterative) TGFM tissue-specific prior')
variant_prob_bootstrapped_distr, tissue_probs_bootstrapped_distr = learn_iterative_variant_gene_tissue_prior_pip_level_bootstrapped(component_level_summary_file, tissue_names, per_window_abf_output_stem, max_iter=args.num_iterations, n_bootstraps=args.num_posterior_samples)


# Print to output
variant_gene_distr_prior_output_file = args.out + '_tissue_specific_prior_summary.txt'
t = open(variant_gene_distr_prior_output_file,'w')
t.write('element_name\tmean_prior\tbootstrapped_prior_distribution\n')
t.write('variant\t' + str(np.mean(variant_prob_bootstrapped_distr)) + '\t' + ';'.join(variant_prob_bootstrapped_distr.astype(str)) + '\n')
for tiss_iter, tissue_name in enumerate(tissue_names):
	t.write(tissue_name + '\t' + str(np.mean(tissue_probs_bootstrapped_distr[tiss_iter,:])) + '\t' + ';'.join(tissue_probs_bootstrapped_distr[tiss_iter,:].astype(str)) + '\n')
t.close()




# Delete unneccessary files
os.system('rm ' + component_level_summary_file)
os.system('rm ' + concatenated_tgfm_no_sampling_pip_summary_file)
os.system('rm ' + per_window_abf_output_stem + '*')

print('DONE')





