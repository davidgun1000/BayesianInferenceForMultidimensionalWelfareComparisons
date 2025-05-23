This folder contains Matlab function to estimate four dimensional well-being 
attributes using Gaussian copula models. See Gunawan, Griffiths, and Chotikapanich
for more details. The users need to supply their own datasets.

The main function is main_prog_estimation.m
This function is applied to estimate the multivariate well-being distributions for the years
2001, 2010, 2015, and 2019 separately. 

param2001.mat contains parameter estimates for the year 2001.
param2010.mat contains parameter estimates for the year 2010.
param2015.mat contains parameter estimates for the year 2015.
param2019.mat contains parameter estimates for the year 2019.

List of functions:
betaden.m: compute beta pdf.
betalogpdf.m: compute the log density of beta distribution.
betarandom.m: generate random numbers from beta distribution.
dirich_rnd.m: generate random samples from Dirichlet distribution.
gamm_rnd.m: generate random samples from gamma distribution, required by dirich_rnd.m
generate_pseudo_representative.m: generate pseudo representative sample using Bayesian bootstrap procedure.
IG_PDF_used.m: compute the PDF of inverted Gamma distribution.
InitialSamplerCopula.m: function to generate the initial values of the model parameters.
log_IG_PDF_used.m: compute the log of inverted gamma density.
logdir3pdf_general.m: compute the log of dirichlet density.
MainSamplerCopula.m: generate model parameters using MCMC.
sample_categorical_parameter.m: generate MCMC samples for the education and happiness distribution.
sampling_beta_parameter.m: generate MCMC samples for the mental health distribution, assumed to follow mixture of beta distributions.
sampling_gamma_parameter.m: generate MCMC samples for the income distribution, assuming to follow mixture of gamma densities. 
sampling_copula_parameter.m: generate MCMC samples for the correlation parameters in the Gaussian copula.