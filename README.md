# Bayesian Adaptive Scalar-on-function Regression (BASOFR)

# Author Contributions Checklist Form


## Data

### Abstract 
This research uses a Bayesian scalar-on-function regression (SOFR) model to analyze the association between prenatal air pollution exposure and later educational outcomes of a large cohort of North Carolina (NC) students. To model the high-resolution and high-dimensional air pollution data, we develop a new Bayesian adaptive scalar-on-function regression (BASOFR) model that features dynamic shrinkage processes, scalable computing, and novel decision analysis tools that allow us to extract the critical time periods of prenatal exposure – if any – that are predictive of adverse educational outcomes.

### Availability 

The NC data cannot be released due to privacy protections; however, access to the data can occur through establishing affiliation with the Children’s Environmental Health Initiative. 

### Simulated data
The code for producing all simulated data and simulation results is included here. 

## Code

This folder contains all necessary files to fit the BASOFR model as well as reproduce all the simulated data and simulation results.  The contents of this folder are as follows:

* ./main_R_code/fit_BASOFR_Sim.R: code for producing the simulated data and fitting the BASOFR model with it
* ./main_R_code/Replicate_Figure3.R: code for reproducing the result shown in Figure 3 of the original paper
* ./main_R_code/Replicate_SimulationResults.R: code for reproducing the simulation results shown in Section 5 and supplementary material of the original paper
* ./main_R_code/MCMC_BASOFR_n500_SNR5_smooth.RDATA: MCMC chain to reproduce Figure 3
* ./helper_functions/dhs_sampling.R: functions that update the parameters of the dynamic horseshoe shrinkage process(DHS) 
* ./helpter_functions/Sim_BASOFR_BLISS_LP_P.R: function that runs simulations using BASOFR and all competitors shown in the original paper
* README.md: README file containing instructions on how to run the code

fit_BASOFR_Sim.R is the main code for fitting the model.  This code uses the following libraries:

* fda 5.4.0
* BayesLogit 2.1
* ggplot2 3.3.5
* latex2exp 0.5.0
* MASS 7.3-54

To run the code, the user will need to specify the following variables (all specified at the top of the file):

* n - the sample size
* STN - the signal-to-noise ratio 
* smoothness - the smoothness of the functional predictors

To run a MCMC with 20000 iterations, fit_BASOFR_Sim.R takes about 40 seconds for sample size 50, 40 seconds for sample size 100, 50 seconds for sample size 500, 

on a MacBook Pro, 2.6GHz Intel Core i7

