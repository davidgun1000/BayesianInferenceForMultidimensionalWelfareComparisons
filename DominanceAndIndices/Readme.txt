This folder contains code to compute Multidimensional indices and dominance probabilities. 
The main function is main_prog_dominance.m for computing indices and dominance probabilities.

param2001.mat contains the model parameters for the year 2001.
param2010.mat contains the model parameters for the year 2010.
param2015.mat contains the model parameters for the year 2015.
param2019.mat contains the model parameters for the year 2019.

The main function calls a pairs of files above, for example param2001.mat and
param2010.mat, to compute the multidimensional indices for the years 2001 and 2010 and
to obtain dominance probabilities between years 2001 and 2010. 

List of functions:
calc_corr.m: compute Pearson and Kendall tau correlations. 
calc_mean_median_educhappiness.m: compute mean and median values of education and happiness distributions.
calc_MWI.m: calculate Multivariate welfare indices for each year.
computing_Cowell_idx.m: compute cowell index for education and happiness distributions.
generate_pseudo_representative: to generate pseudo representative sample using Bayesian bootstrap method of Gunawan et  al 2020
main_prog_dominance.m: the main function to compute various dominance probabilities and indices
main_prog_dominance_first_v3.m: the main function to compute dominance probabilities related to U1 dominance conditions.