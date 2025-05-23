This folder contains the code to estimate univariate mixture of gamma and beta distributions.
We perform validation set approach with Mean Absolute Error as the loss function to select
the optimal number of components for the mixture of gamma distributions and the 
mixture of beta distributions.

The main function is train_test_data.m

List of functions:
betaden.m: compute beta pdf.
betalogpdf.m: compute the log density of beta distribution.
betarandom.m: generate random numbers from beta distribution.
data_adj_train_test.m: split the data to training and test sets.
dirich_rnd.m: generate random samples from Dirichlet distribution.
gamm_rnd.m: generate random samples from gamma distribution, required by dirich_rnd.m
generate_pseudo_representative.m: generate pseudo representative sample using Bayesian bootstrap procedure.
IG_PDF_used.m: compute the PDF of inverted Gamma distribution.
InitialSamplerBetaMix.m: function to generate the initial values of the mixture of beta parameters.
InitialSamplerGammaMix.m: function to generate the initial values of the mixture of gamma parameters.
log_IG_PDF_used.m: compute the log of inverted gamma density.
logdir3pdf_general.m: compute the log of dirichlet density.
MainSamplerBetaMix.m: generate mixture of beta parameters using MCMC.
MainSamplerGammaMix.m: generate mixture of gamma parameters using MCMC.
sampling_beta_parameter_only_Mix.m: MCMC steps for mixture of beta parameters
sampling_gamma_parameter_only_Mix.m: MCMC steps for the mixture of gamma parameters.

The results are given in the folder "Results".