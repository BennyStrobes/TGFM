import sys
import pandas as pd
import numpy as np 
import os 
import pdb
import scipy.special
import time



class TGFM(object):
	def __init__(self, L=10, gene_init_log_pi=None, variant_init_log_pi=None):
		# Prior on gamma distributions defining residual variance and
		self.L = L
		self.gene_init_log_pi = gene_init_log_pi
		self.variant_init_log_pi = variant_init_log_pi
	def fit(self, twas_data_obj):
		""" Fit the model.
			Args:
			twas_data_obj
		"""

		print('###############################')
		print('TGFM without sampling')
		print('###############################')

		#####################
		# Initialize variables
		self.initialize_variables(twas_data_obj)

		return

	def initialize_variables(self, twas_data_obj):
		# Number of genes
		self.G = len(twas_data_obj['gene_tissue_pairs'])
		# Number of variants
		self.K = len(twas_data_obj['variants'])
		# Gene names
		self.gene_tissue_pairs = twas_data_obj['gene_tissue_pairs']
		# GWAS sample size
		self.NN = twas_data_obj['gwas_sample_size']

		
		# Compute nominal twas z-scores
		self.nominal_twas_z = np.zeros(self.G)


		# Initialize variational distributions defining alphas (the causal effect of genetically-predicted expression in each gene on the trait)
		# Currently using null intitialization
		self.alpha_mu = np.zeros((self.L, self.G))
		self.alpha_var = np.ones((self.L, self.G))
		self.alpha_phi = np.ones((self.L, self.G))/(self.G + self.K)
		self.component_variances = np.ones(self.L)*1e4


		self.prev_alpha_mu = np.copy(self.alpha_mu)

		# Initialize variational distribution defining betas (the causal effect of pleiotropic genotype on the trait)
		self.beta_mu = np.zeros((self.L, self.K))
		self.beta_var = np.ones((self.L, self.K))
		self.beta_phi = np.ones((self.L, self.K))/(self.G + self.K)


		# Initialize KL terms
		self.KL_terms = np.zeros(self.L)
		self.LBF_terms = np.zeros(self.L)
		self.elbo = 0.0

		return
