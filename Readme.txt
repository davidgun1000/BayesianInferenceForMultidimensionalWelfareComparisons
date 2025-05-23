This MATLAB software implements the MCMC sampler for estimating the copula models for modeling four well-being attributes: income, mental health, education,
and happiness variables and the algorithm to compute different indices and probabilities of dominance. 
from the paper "Bayesian Inference from Multidimensional Welfare Comparisons".
There are four folders:
1) ChoosingNumberOfComponents: This folder contains code to choose number of components of gamma and beta mixture models for income and mental health, respectively.
2) Estimation: This folder contains code to implement MCMC samplers for estimating copula models for modeling joint distribution of well-being attributes. 
3) DominanceandIndices: This folder contains code to estimate different indices and dominance probabilities. 
4) Simulation: This folder contains code to estimate the four dimensional copula models using simulated data with equally weighted observations.

Folder 2001_2010 contains the dominance probabilities and indices for comparing 2001 and 2010.
Folder 2010_2015 contains the dominance probabilities and indices for comparing 2010 and 2015
Folder 2015_2019 contains the dominance probabilities and indices for comparing 2015 and 2019

The real data were obtained from four waves of the Household, Income and Labour Dynamics in Australia (HILDA) survey (Release 19) corresponding to years 2001,
2010, 2015, and 2019. Prospective users can request the data at https://melbourneinstitute.unimelb.edu.au/hilda/for-data-users#accessing and use the data after satisfying the confidentiality requirements.

The HILDA survey is a national representative longitudinal survey which began in Australia in 2001
and is conducted annually. It was initiated and is funded by the Australian Government Department of
Social Services (DSS) and is managed by the Melbourne Institute of Applied Economics and Social
Research (Melbourne Institute). Full details are provided on its website which has now been referenced
in the paper (https://melbourneinstitute.unimelb.edu.au/hilda). Permission is required to use the data,
but ethics approval is not. Data on key variables concerning family and household structure, as well as
data on education, income, health, life satisfaction and other variables relating to economic and
subjective wellbeing are collected. It is a longitudinal sample, but individuals do leave and enter the
panel. 
