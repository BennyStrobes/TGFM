import sys
import numpy as np 
import os 
import pdb
import scipy.special



def fill_in_causal_effect_size_matrix(null_mat, bs_eqtls_pmces_sparse):
	null_mat[bs_eqtls_pmces_sparse[:,0].astype(int), bs_eqtls_pmces_sparse[:,1].astype(int)] = bs_eqtls_pmces_sparse[:,2]
	return null_mat


def compute_expected_pips_from_sampler_phis(gene_phis):
	# Number of samples
	n_bs = gene_phis[0].shape[0]
	# Number of genes
	n_genes = gene_phis[0].shape[1]
	# Number of components
	LL = len(gene_phis)

	# Initialize per-sample gene PIPs
	gene_pips = np.ones((n_bs, n_genes))
	# Fill in per sample gene PIPs
	for component_iter in range(LL):
		gene_pips = gene_pips*(1.0 - gene_phis[component_iter])
	gene_pips = 1.0 - gene_pips

	# Compute Expected value across samples
	expected_gene_pips = np.mean(gene_pips,axis=0)

	return expected_gene_pips


class TGFM(object):
	def __init__(self, L=10, max_iter=5, gene_log_pi=None, variant_log_pi=None):
		self.L = L  # Number of components
		self.max_iter = max_iter  # Total iterations
		self.gene_log_pi = gene_log_pi  # Log prior probabilities for gene-tissue pairs
		self.variant_log_pi = variant_log_pi  # Log prior probabilities for non-mediated variants
	def fit(self, tgfm_data_obj, phi_init=None, mu_init=None, mu_var_init=None):
		""" Fit the model.
			Args:
			tgfm_data_obj
		"""

		#####################
		# Initialize variables
		print('Initialize variables')
		self.initialize_variables(tgfm_data_obj, phi_init, mu_init, mu_var_init)

		print('Iterative optimization')
		# Begin iterative inference
		for itera in range(self.max_iter):
			# Update alpha (ie predicted effect of each gene on the trait) and beta (ie predicted effect of variant on trait)
			self.update_susie_effects(self.gene_log_pi, self.variant_log_pi, self.component_variances, self.component_variances)

			# Update variance of each component in each bootstrap/posterior sample
			self.update_component_variances()

			# Keep track of the number of iterations
			self.iter = self.iter + 1

		# COMPUTE PIPS for each bootstrap
		print('Compute pips')
		self.compute_pips()  # Compute PIPs for gene-tissue pairs and non-mediated variants
		self.compute_gene_pips()  # Compute Gene PIPs

		return

	def compute_gene_pips(self):
		# Extract unique gene names
		unique_gene_names = np.unique(self.gene_names)

		# Create mapping from gene name to all gene-tissue pair level indices that correspond to the gne
		gene_to_indices = {}
		for ii, gene_name in enumerate(self.gene_names):
			if gene_name not in gene_to_indices:
				gene_to_indices[gene_name] = []
			gene_to_indices[gene_name].append(ii)

		# Total number of unique genes
		n_genes = len(gene_to_indices)
		# Ordered list of genes
		ordered_genes = [*gene_to_indices]
		# Number of components
		LL = len(self.alpha_phis)

		# Array to keep track of gene_phis (gene inclusion probabilities)
		gene_phis = []

		# Loop through components
		for ll in range(LL):
			# Initialize gene_phi for this component
			gene_phi = np.zeros((self.alpha_phis[ll].shape[0], n_genes))

			# Loop through unique genes
			for ii,gene_name in enumerate(ordered_genes):
				gene_phi[:, ii] = np.sum(self.alpha_phis[ll][:, gene_to_indices[gene_name]],axis=1)
			gene_phis.append(gene_phi)

		# Compute per-sample gene PIPs, and then average across samples
		expected_gene_pips= compute_expected_pips_from_sampler_phis(gene_phis)

		# Create mapping from gene name to gene
		self.gene_name_to_gene_pip = {}
		for ii, gene_name in enumerate(ordered_genes):
			self.gene_name_to_gene_pip[gene_name] = expected_gene_pips[ii]

		return


	def compute_pips(self):
		self.alpha_pips = np.ones((self.n_bs, self.G))
		self.beta_pips = np.ones((self.n_bs, self.K))

		for component_iter in range(self.L):
			self.alpha_pips = self.alpha_pips*(1.0 - self.alpha_phis[component_iter])
			self.beta_pips = self.beta_pips*(1.0 - self.beta_phis[component_iter])

		self.alpha_pips = 1.0 - self.alpha_pips
		self.beta_pips = 1.0 - self.beta_pips

		self.expected_alpha_pips = np.mean(self.alpha_pips,axis=0)
		self.expected_beta_pips = np.mean(self.beta_pips,axis=0)

		return




	def update_susie_effects(self, expected_log_pi, expected_log_variant_pi, alpha_component_variances, beta_component_variances):
		
		# Initialize array to keep track of causal eqtl effect sizes across bootstraps
		gene_eqtl_pmces = np.zeros((self.G, self.K))
		
		# Loop through components
		for l_index in range(self.L):
			if l_index == 0:
				prev_l_index = self.L - 1
			else:
				prev_l_index = l_index - 1

			# Include current component (as it is currently removed)
			component_non_med_pred = self.beta_mus[l_index]*self.beta_phis[l_index] - self.beta_mus[prev_l_index]*self.beta_phis[prev_l_index]  # Predicted non-mediated effects
			self.residual = self.residual + np.dot(self.gene_trait_pred[l_index] - self.gene_trait_pred[prev_l_index] + component_non_med_pred, self.srs_inv)  # Residual with both variant-trait effects and variant-gene-trait effects removed
			self.residual_incl_alpha = self.residual_incl_alpha + np.dot(component_non_med_pred, self.srs_inv)  # Residual with only variant-trait effects removed
			
			# Loop through samples/bootstraps
			for bs_iter in range(self.n_bs):

				# Load in gene eQTL effects for this bootstrap
				gene_eqtl_pmces = gene_eqtl_pmces*0.0
				bs_eqtls_pmces_sparse = self.sparse_sampled_gene_eqtl_pmces[bs_iter]
				gene_eqtl_pmces = fill_in_causal_effect_size_matrix(gene_eqtl_pmces, bs_eqtls_pmces_sparse)

				#################
				# Alpha effects  (gene trait effects)
				#################
				# Calculate terms inolved in update
				gene_weights = self.agg_gene_trait[bs_iter] - (self.alpha_mus[l_index][bs_iter,:]*self.alpha_phis[l_index][bs_iter,:])
				b_terms = np.dot(gene_eqtl_pmces, np.multiply(self.residual_incl_alpha[bs_iter,:], self.s_inv_2_diag)) + np.dot(self.precomputed_gene_gene_terms[bs_iter], gene_weights)
				a_terms = self.precomputed_a_terms[bs_iter] - .5*(1.0/alpha_component_variances[bs_iter][l_index])
				# Perform update
				mixture_alpha_var = -1.0/(2.0*a_terms)
				mixture_alpha_mu = b_terms*mixture_alpha_var

				#################
				# Beta effects (variant trait effects)
				#################
				# Calculate terms involved in update
				variant_b_terms = self.residual[bs_iter,:]*self.s_inv_2_diag
				variant_a_terms = (-.5*self.D_diag) - .5*(1.0/beta_component_variances[bs_iter][l_index])
				# Perform update
				mixture_beta_var = -1.0/(2.0*variant_a_terms)
				mixture_beta_mu = variant_b_terms*mixture_beta_var

			
				################
				# Normalization (across beta and alpha)
				###############
				un_normalized_lv_alpha_weights = expected_log_pi[:,bs_iter] - (.5*np.log(alpha_component_variances[bs_iter][l_index])) + (.5*np.square(mixture_alpha_mu)/mixture_alpha_var) + (.5*np.log(mixture_alpha_var))
				un_normalized_lv_beta_weights = expected_log_variant_pi[:,bs_iter] - (.5*np.log(beta_component_variances[bs_iter][l_index])) + (.5*np.square(mixture_beta_mu)/mixture_beta_var) + (.5*np.log(mixture_beta_var))
				normalizing_term = scipy.special.logsumexp(np.hstack((un_normalized_lv_alpha_weights, un_normalized_lv_beta_weights)))


				# Update agg_gene_trait pt 1
				self.agg_gene_trait[bs_iter] = self.agg_gene_trait[bs_iter] - (self.alpha_mus[l_index][bs_iter,:]*self.alpha_phis[l_index][bs_iter,:])

				# Save results to global model parameters
				mixture_alpha_phi = np.exp(un_normalized_lv_alpha_weights-normalizing_term)
				self.alpha_phis[l_index][bs_iter,:] = mixture_alpha_phi
				self.alpha_mus[l_index][bs_iter,:] = mixture_alpha_mu
				self.alpha_vars[l_index][bs_iter,:] = mixture_alpha_var

				mixture_beta_phi = np.exp(un_normalized_lv_beta_weights-normalizing_term)
				self.beta_phis[l_index][bs_iter,:] = mixture_beta_phi
				self.beta_mus[l_index][bs_iter,:] = mixture_beta_mu
				self.beta_vars[l_index][bs_iter,:] = mixture_beta_var


				# Update KL Terms
				lbf = np.hstack((un_normalized_lv_alpha_weights-expected_log_pi[:, bs_iter], un_normalized_lv_beta_weights-expected_log_variant_pi[:, bs_iter]))
				# Save LBF
				self.alpha_lbfs[l_index][bs_iter,:] = un_normalized_lv_alpha_weights-expected_log_pi[:, bs_iter]
				self.beta_lbfs[l_index][bs_iter,:] = un_normalized_lv_beta_weights-expected_log_variant_pi[:, bs_iter]
				# Update KL terms
				maxlbf = np.max(lbf)
				ww = np.exp(lbf - maxlbf)
				ww_weighted = ww*np.exp(np.hstack((expected_log_pi[:, bs_iter], expected_log_variant_pi[:, bs_iter])))
				kl_term1 = -(np.log(np.sum(ww_weighted)) + maxlbf) # THIS TERM IS CORRECT
				kl_term2 = np.sum((mixture_beta_mu*mixture_beta_phi)*variant_b_terms) + np.sum((mixture_alpha_mu*mixture_alpha_phi)*b_terms)
				kl_term3 = -.5*(np.sum((self.NN - 1)*mixture_beta_phi*(np.square(mixture_beta_mu) + mixture_beta_var)) + np.sum((self.NN - 1)*mixture_alpha_phi*(np.square(mixture_alpha_mu) + mixture_alpha_var)))
				self.KL_terms[l_index][bs_iter] = kl_term1 + kl_term2 + kl_term3

				# Update agg_gene_trait pt 2
				self.agg_gene_trait[bs_iter] = self.agg_gene_trait[bs_iter] + (self.alpha_mus[l_index][bs_iter,:]*self.alpha_phis[l_index][bs_iter,:])

				# Update gene trait pred
				self.gene_trait_pred[l_index][bs_iter,:] = np.dot((self.alpha_mus[l_index][bs_iter,:])*(self.alpha_phis[l_index][bs_iter,:]), gene_eqtl_pmces)
		return


	def update_component_variances(self):
		# NOTE: COMPONENT VARIANCE IS SHARED ACROSS BETA AND ALPHA: Maybe a bad idea??
		for bs_iter in range(self.n_bs):
			for l_iter in range(self.L):
				self.component_variances[bs_iter][l_iter] = np.sum((np.square(self.alpha_mus[l_iter][bs_iter,:]) + self.alpha_vars[l_iter][bs_iter,:])*self.alpha_phis[l_iter][bs_iter,:]) + np.sum((np.square(self.beta_mus[l_iter][bs_iter,:]) + self.beta_vars[l_iter][bs_iter,:])*self.beta_phis[l_iter][bs_iter,:])

		return



	def initialize_variables(self, tgfm_data_obj, phi_init, mu_init, mu_var_init):
		# Number of genes
		self.G = len(tgfm_data_obj['gene_tissue_pairs'])
		# Number of variants
		self.K = len(tgfm_data_obj['variants'])
		# Gene-tissue names names
		self.gene_tissue_pairs = tgfm_data_obj['gene_tissue_pairs']
		# Gene names
		self.gene_names = tgfm_data_obj['genes']
		# GWAS sample size
		self.NN = tgfm_data_obj['gwas_sample_size']

		# Get gwas z-scores
		self.gwas_variant_z = tgfm_data_obj['gwas_beta']/tgfm_data_obj['gwas_beta_se']

		# Generate S matrix (Diagonail matrix of gwas variances)
		s_squared_vec = np.square(tgfm_data_obj['gwas_beta_se']) + (np.square(tgfm_data_obj['gwas_beta'])/tgfm_data_obj['gwas_sample_size'])
		s_vec = np.sqrt(s_squared_vec)
		S_mat = np.diag(s_vec)
		S_inv_mat = np.diag(1.0/s_vec)
		S_inv_2_mat = np.diag(1.0/np.square(s_vec))

		# Compute (S^-1)R(S^-1) taking advantage of fact that S^-1 is a diagonal matrix
		D_mat = np.multiply(np.multiply(np.diag(S_inv_mat)[:, None], tgfm_data_obj['reference_ld']), np.diag(S_inv_mat))
		# Compute (S)R(S^-1) taking advantage of fact that S and S^-1 is a diagonal matrix
		srs_inv_mat = np.multiply(np.multiply(np.diag(S_mat)[:, None], tgfm_data_obj['reference_ld']), np.diag(S_inv_mat))

		# Generate data object containing statistics that are precomputed
		self.srs_inv = srs_inv_mat
		self.s_inv_2_diag = np.diag(S_inv_2_mat)
		self.D_diag = np.diag(D_mat)

		# Number of posterior samples/bootstrapped
		self.n_bs = len(tgfm_data_obj['sparse_sampled_gene_eqtl_pmces'])

		# Precompute some terms related to gene variances and gene-gene covariance that
		self.precomputed_a_terms = []  # Related to gene variance
		self.precomputed_gene_gene_terms = []  # Related to gene-gene covariance
		self.agg_gene_trait = []
		# Also compute nominal twas z scores
		self.nominal_twas_z = []
		# Also initialize component variance parameters
		self.component_variances = []

		# Initialize matrix of eqtl effects (will be filled in for each sample/bootstrap in loop below)
		bs_eqtls_pmces = np.zeros((self.G, self.K))

		# Loop through samples/bootstraps
		for bs_iter in range(self.n_bs):
			# Get Causal eqtl effect sizes for this sample/bootstrap
			bs_eqtls_pmces = bs_eqtls_pmces*0.0
			bs_eqtls_pmces_sparse = tgfm_data_obj['sparse_sampled_gene_eqtl_pmces'][bs_iter]
			bs_eqtls_pmces = fill_in_causal_effect_size_matrix(bs_eqtls_pmces, bs_eqtls_pmces_sparse)

			# Precompute terms related gene variances and gene-gene covariance in this bootstrap
			# Only perform computation for snps with non-zero eqtl effect in at least one gene (sppeds it up a bit particular cause eqtl effects are very spare)
			non_zero_snp_indices = np.sum(bs_eqtls_pmces!=0.0,axis=0) != 0.0
			bs_precomputed_gene_gene_terms = -np.dot(np.dot(bs_eqtls_pmces[:, non_zero_snp_indices],D_mat[:, non_zero_snp_indices][non_zero_snp_indices,:]), np.transpose(bs_eqtls_pmces[:, non_zero_snp_indices]))
			bs_precomputed_a_terms = .5*np.diag(bs_precomputed_gene_gene_terms)
			self.precomputed_a_terms.append(bs_precomputed_a_terms)
			self.precomputed_gene_gene_terms.append(bs_precomputed_gene_gene_terms)

			# Initialize variance of component in this bootstrap
			self.component_variances.append(np.ones(self.L)*1e4)

			# Initialize aggregrate gene trait effect
			self.agg_gene_trait.append(np.zeros(self.G))

			# Compute nominal twas zs
			bs_nominal_twas_z = np.dot(bs_eqtls_pmces, self.gwas_variant_z)
			self.nominal_twas_z.append(bs_nominal_twas_z)

		# Initialize all other model parameters
		self.alpha_mus = []
		self.alpha_vars = []
		self.alpha_phis = []
		self.beta_mus = []
		self.beta_vars = []
		self.beta_phis = []
		self.gene_trait_pred = []
		self.KL_terms = []
		self.alpha_lbfs = []
		self.beta_lbfs = []

		if phi_init is not None and mu_init is not None and mu_var_init is not None:
			# Pre-specified initialization
			# Set parameters in all bootstraps to phi_init, mu_init, and mu_var_init
			for l_iter in range(self.L):
				# Initialize variational distributions defining alphas (the causal effect of genetically-predicted expression in each gene on the trait)
				self.alpha_mus.append(np.repeat(np.reshape(mu_init[l_iter,:],(1, mu_init.shape[1])),self.n_bs,axis=0)[:,:self.G])
				self.alpha_vars.append(np.repeat(np.reshape(mu_var_init[l_iter,:],(1, mu_var_init.shape[1])),self.n_bs,axis=0)[:,:self.G])
				self.alpha_phis.append(np.repeat(np.reshape(phi_init[l_iter,:],(1, phi_init.shape[1])),self.n_bs,axis=0)[:,:self.G])

				# Initialize variational distribution defining betas (the causal effect of pleiotropic genotype on the trait)
				self.beta_mus.append(np.repeat(np.reshape(mu_init[l_iter,:],(1, mu_init.shape[1])),self.n_bs,axis=0)[:,self.G:])
				self.beta_vars.append(np.repeat(np.reshape(mu_var_init[l_iter,:],(1, mu_var_init.shape[1])),self.n_bs,axis=0)[:,self.G:])
				self.beta_phis.append(np.repeat(np.reshape(phi_init[l_iter,:],(1, phi_init.shape[1])),self.n_bs,axis=0)[:,self.G:])
				
				# Initialize quantity to keep track of predicted gene-trait effects
				self.gene_trait_pred.append(np.zeros((self.n_bs, self.K)))
				if np.sum(np.sum(np.repeat(np.reshape(mu_init[l_iter,:],(1, mu_init.shape[1])),self.n_bs,axis=0)[:,:self.G])) != 0.0:
					print('assumption erroror: non-zedro initialized gene effect')
					pdb.set_trace()

				# Initialize KL terms
				self.KL_terms.append(np.zeros(self.n_bs))
				# Initialize lbf
				self.alpha_lbfs.append(np.zeros((self.n_bs, self.G)))
				self.beta_lbfs.append(np.zeros((self.n_bs, self.K)))
		else:
			# Null initialization
			for l_iter in range(self.L):
				# Initialize variational distributions defining alphas (the causal effect of genetically-predicted expression in each gene on the trait)
				# Currently using null intitialization
				self.alpha_mus.append(np.zeros((self.n_bs, self.G)))
				self.alpha_vars.append(np.ones((self.n_bs, self.G)))
				self.alpha_phis.append(np.ones((self.n_bs, self.G))/(self.G + self.K))

				# Initialize variational distribution defining betas (the causal effect of pleiotropic genotype on the trait)
				self.beta_mus.append(np.zeros((self.n_bs, self.K)))
				self.beta_vars.append(np.ones((self.n_bs, self.K)))
				self.beta_phis.append(np.ones((self.n_bs, self.K))/(self.G + self.K))

				# Initialize quantity to keep track of predicted gene-trait effects
				self.gene_trait_pred.append(np.zeros((self.n_bs, self.K)))

				# Initialize KL terms
				self.KL_terms.append(np.zeros(self.n_bs))
				# Initialize lbf
				self.alpha_lbfs.append(np.zeros((self.n_bs, self.G)))
				self.beta_lbfs.append(np.zeros((self.n_bs, self.K)))

		# Add sparse eqtl effects sizes to data object
		self.sparse_sampled_gene_eqtl_pmces = tgfm_data_obj['sparse_sampled_gene_eqtl_pmces']

		# Keep track of iterations
		self.iter = 0

		# Compute residual gwas effect
		# For computation purposes: one with gene-trait effects (alphas) and one without
		self.residual = []
		self.residual_incl_alpha = []
		for bs_iter in range(self.n_bs):
			self.residual.append(tgfm_data_obj['gwas_beta'])
			self.residual_incl_alpha.append(tgfm_data_obj['gwas_beta'])
		self.residual = np.asarray(self.residual)
		self.residual_incl_alpha = np.asarray(self.residual_incl_alpha)

		# Remove initialized non-mediated effects
		# Currently not implemented to remove initialized mediated efffects (but will return an  errror if we try to)
		for l_index in range(self.L - 1):
			component_non_med_pred = self.beta_mus[l_index]*self.beta_phis[l_index]
			self.residual = self.residual - np.dot(component_non_med_pred, self.srs_inv)
			self.residual_incl_alpha = self.residual_incl_alpha - np.dot(component_non_med_pred, self.srs_inv)



		return




