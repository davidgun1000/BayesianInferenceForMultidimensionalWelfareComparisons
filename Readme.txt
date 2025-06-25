This MATLAB package implements a Markov Chain Monte Carlo (MCMC) sampler for estimating copula-based models used to jointly model four well-being attributes: income, mental health, education, and happiness. It also includes algorithms for computing various inequality and welfare indices as well as posterior probabilities of Lorenz and stochastic dominance.

The methodology and models are described in the paper:
"Bayesian Inference from Multidimensional Welfare Comparisons"
by David Gunawan, William Griffiths, and Duangkamon Chotikapanich.

Folder Structure

ChoosingNumberOfComponents/
Contains MATLAB code for selecting the number of components in the mixture models:
Gamma mixture model for income
Beta mixture model for mental health

Estimation/
Includes MCMC routines for estimating the copula models that capture the joint distribution of the four well-being attributes.

DominanceandIndices/
Provides code to compute:
Various welfare and inequality indices
Posterior probabilities of dominance between different time periods or subgroups

Simulation/
Contains code for estimating the four-dimensional copula model using simulated data, where all observations are equally weighted.

2001_2010/, 2010_2015/, 2015_2019/
These folders include the output of the analysis:
Estimated indices
Dominance probabilities for comparing the years 2001 vs. 2010, 2010 vs. 2015, and 2015 vs. 2019

Data Source
The empirical analysis is based on data from four waves of the Household, Income and Labour Dynamics in Australia (HILDA) Survey, Release 19, corresponding to the years 2001, 2010, 2015, and 2019.
The HILDA Survey is a nationally representative longitudinal household study that began in 2001. It collects detailed information on:
Family and household structure 
Education and income
Health and mental well-being
Life satisfaction and other aspects of subjective and economic well-being
The survey is:
Funded by the Australian Government Department of Social Services (DSS)
Managed by the Melbourne Institute of Applied Economic and Social Research
More information is available at: https://melbourneinstitute.unimelb.edu.au/hilda
Prospective users must request access to the data and meet confidentiality requirements. However, ethics approval is not required for use.
