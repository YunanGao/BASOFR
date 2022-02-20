#############################################
## Code to Replicate the simulation result ##
#############################################
## Load Libraries
library(bliss)
library(BayesLogit)

## Source in Helper Functions
source("../helper_functions/dhs_sampling.R")
source("../helper_functions/Sim_BASOFR_BLISS_LP_P.R")

## All the simulation designs we considered
n.array = c(50,50,100,100,500,500,10000,10000)
STN.array = c(10,10,7,7,5,5,0.5,0.5)
smoothness.array = rep(c("smooth","nonsmooth"),4)
simulation.setting = as.data.frame(cbind(n.array, STN.array, smoothness.array))
colnames(simulation.setting) = c("n","STN","smoothness")
simulation.setting$n = as.numeric(simulation.setting$n)
simulation.setting$STN = as.numeric(simulation.setting$STN)
simulation.setting

## We compared the BASOFR with the local P-spline, P-spline method, and the BLISS method (if sample size <= 100)
## in terms of L2 error, mean credible interval width and empirical pointwise coverage for the 95% CIs computed 
## under each model. The results are shown in the paper. 

## To reproduce the result, input the sample size, signal-to-noise ratio, and smoothness

## for example
n= 50; SNR = 10; mode = "smooth"

# The argument N is the number of simulations, in the original paper N=50
test = Sim_BASOFR_LP_P_BLISS(N=50, n=n, STN = SNR, smoothness = mode)

result = as.data.frame(test)

## L2 error boxplot
par(mai = c(1, 2, 1, 0.5), las = 1);
if(n <= 100){
  L2_beta = cbind(result$L2_BASOFR,result$L2_LP, result$L2_P, result$L2_BLISS)
  bp0 = boxplot(L2_beta, plot = F)
  bp = boxplot(L2_beta, 
               main = paste('L2 error:\n  n =', n, ', SNR =', SNR, ',', mode  ),
               ylim = range(1.15*bp0$stats[5,], .85*bp0$stats[2,], na.rm=TRUE),
               names=c("BASOFR", "Local P-spline", "P-spline", "Bliss"),
               horizontal = TRUE,
               col = 'gray',
               cex.lab = 2, cex.main = 2, cex.axis = 1.5, notch = TRUE, outline = FALSE);
  
}else{
  L2_beta = cbind(result$L2_BASOFR,result$L2_LP, result$L2_P)
  bp0 = boxplot(L2_beta, plot = F)
  bp = boxplot(L2_beta, 
               main = paste('L2 error:\n  n =', n, ', SNR =', SNR, ',', mode  ),
               ylim = range(1.15*bp0$stats[5,], .85*bp0$stats[2,], na.rm=TRUE),
               names=c("BASOFR", "Local P-spline", "P-spline"),
               horizontal = TRUE,
               col = 'gray',
               cex.lab = 2, cex.main = 2, cex.axis = 1.5, notch = TRUE, outline = FALSE);
}


## mean interval widths boxplot
if(n<=100){
  ci_width_beta = cbind(result$area_BASOFR, result$area_LP, result$area_P, result$area_BLISS)
  ci_cover_beta = cbind(result$cov_BASOFR, result$cov_LP, result$cov_P, result$cov_BLISS)
  ci_cover_beta_median  = apply(ci_cover_beta, 2, median)
  bp0 = boxplot(ci_width_beta, plot = F)
  #n=100;SNR=7
  bp = boxplot(ci_width_beta, 
               main = paste('Mean interval widths:\n  n =', n,  ', SNR =',  SNR, ',',  mode),
               ylim = range(1.20*bp0$stats[5,], .85*bp0$stats[2,], na.rm=TRUE),
               names=c("BASOFR", "Local P-spline", "P-spline", "Bliss"),
               horizontal = TRUE,
               col = 'gray',
               cex.lab = 2, cex.main = 2, cex.axis = 1.5, notch = TRUE, outline = FALSE);
  text( 
    x=1*bp0$stats[5,],
    y=c(1:4), 
    paste(' ', round(ci_cover_beta_median*100), '%', sep=''),
    cex = 2.5, adj = c(0,NA), col='blue'
  )
}else{
  ci_width_beta = cbind(result$area_BASOFR, result$area_LP, result$area_P)
  ci_cover_beta = cbind(result$cov_BASOFR, result$cov_LP, result$cov_P)
  ci_cover_beta_median  = apply(ci_cover_beta, 2, median)
  bp0 = boxplot(ci_width_beta, plot = F)
  #n=100;SNR=7
  bp = boxplot(ci_width_beta, 
               main = paste('Mean interval widths:\n  n =', n,  ', SNR =',  SNR, ',',  mode),
               ylim = range(1.20*bp0$stats[5,], .85*bp0$stats[2,], na.rm=TRUE),
               names=c("BASOFR", "Local P-spline", "P-spline"),
               horizontal = TRUE,
               col = 'gray',
               cex.lab = 2, cex.main = 2, cex.axis = 1.5, notch = TRUE, outline = FALSE);
  text( 
    x=1*bp0$stats[5,],
    y=c(1:4), 
    paste(' ', round(ci_cover_beta_median*100), '%', sep=''),
    cex = 2, adj = c(0,NA), col='blue'
  )
}


