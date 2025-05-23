This folder contains code for replicating the simulation study described in section I of the online supplement, by assuming equally weighted samples, so we do not need to 
use Bayesian bootstrap method to generate pseudo representative samples. 

The main code is SimulationCopula.m. We have set a seed and the results can be replicated.

The data is simulated data.mat

List of functions:
betaden.m: compute beta pdf.
betalogpdf.m: compute the log density of beta distribution.
betarandom.m: generate random numbers from beta distribution.
dirich_rnd.m: generate random samples from Dirichlet distribution.
gamm_rnd.m: generate random samples from gamma distribution, required by dirich_rnd.m
IG_PDF_used.m: compute the PDF of inverted Gamma distribution.
log_IG_PDF_used.m: compute the log of inverted gamma density.
logdir3pdf_general.m: compute the log of dirichlet density.
sample_categorical_parameter.m: generate MCMC samples for the education and happiness distribution.
sampling_beta_parameter.m: generate MCMC samples for the mental health distribution, assumed to follow mixture of beta distributions.
sampling_gamma_parameter.m: generate MCMC samples for the income distribution, assuming to follow mixture of gamma densities. 
sampling_copula_parameter.m: generate MCMC samples for the correlation parameters in the Gaussian copula.