####################################################
## Code to Fit BASOFR Model to Simulation Data ##
####################################################

###################
## Preliminaries ##
###################

##Load Libraries
library(fda)
library(MASS)
library(BayesLogit)
library(ggplot2)
library(latex2exp)

## Source in Helper Functions
source("../helper_functions/dhs_sampling.R")

## Set sample size, signal-to-noise ratio (SNR)
n=500; STN=5 

## Set grid the functional predictors are observed on
obs.x=101
x = seq(0,1,length.out = obs.x)

## Set true regression function \beta(t)
beta = 8 * (2+exp(20-60*x) + exp(60*x-20))^(-1) - 12*(2 + exp(40-60*x) + exp(60*x-40))^(-1)

#######################
## Simulate the data ##
#######################

## Are the functional predictors smooth or nonsmooth?
smoothness = "smooth"
if(smoothness == "smooth"){l = 0.1}
if(smoothness == "nonsmooth"){l=0.01}

## Simulate the functional observations from a Gaussian process
# set.seed(1) 
x_diff = abs(matrix(rep(x,obs.x), nrow=obs.x)-t(matrix(rep(x,obs.x),nrow=obs.x)))
sigma_x = 0.7
Kernel = sigma_x^2 * exp(-x_diff^2/(2*l^2))
mu = rep(0,obs.x)
x.sim = mvrnorm(n=n, mu=mu, Sigma=Kernel) 

## Simulate the scalar responses 
noise.sd = sd(tcrossprod(x.sim, t(beta))/(obs.x-1))/STN
y = tcrossprod(x.sim, t(beta))/(obs.x-1) + rnorm(n, mean=0,sd=noise.sd)

#######################################################################
## Smooth the functional observations using B-spline basis expansion ##
#######################################################################

## Create the B-spline basis functions
tobs = seq(0,1,0.01);nobs = length(tobs)
knots = c(seq(0,1,0.02));nknots = length(knots)
norder = 4 
nbasis = length(knots) + norder - 2
basis = create.bspline.basis(c(min(tobs), max(tobs)), nbasis = nbasis, 
                             norder=norder, breaks = knots)
basismat = eval.basis(tobs,basis)

## dimension of the B-spline basis (K is also the value of K_B and K_X)
K = ncol(basismat)

## estimate the OLS b-spline basis coefficients
X.star = t(tcrossprod(tcrossprod(solve(crossprod(basismat)), basismat), x.sim))

########################################################
## Transform the SOFR to a multiple linear regression ##
########################################################

## Estimate the J matrix
# (Here, we expand the functional predictors and the regression function on the same B-spline basis )
J = eval.penalty(basisobj= basis, Lfdobj=0, rng=c(0,1))

## the predictors for the multiple linear regression
X.star.star = tcrossprod(X.star, J) 

## design matrix for the multiple linear regression
X = cbind(1, X.star.star) 


###################
## MCMC settings ##
###################
burn <- 10000
iter <- 10000

## total iterations
S = burn + iter 

######################
## Prior parameters ##
######################

## Gamma prior on the SOFR (1/{\sigma^2}) term
a0 = 0.01; b0=0.01 

## Gamma prior on {\lambda_0}^{-2} 
a.lambda=0.01; b.lambda=0.01

#####################
## Starting values ##
#####################

## the intercept term and the B-spline coefficients of the regression function
B=crossprod(t(solve(crossprod(X) + diag(10^(-5), nrow=(K+1)))),crossprod(X, y))

## {lambda_k^2}
D_2 = diff(diag(K), diff=2)
evolParams = initDHS(omega = as.matrix(tcrossprod(D_2, t(as.matrix(B[-1])))))
h_k = as.matrix(evolParams$ht)
lambda_k.sq = as.numeric(exp(h_k)) #lambda_k.sq, for k=2,...,(K_B-1)
lambda0.sq = 1/(rgamma(1,a.lambda, b.lambda)) #lambda_k.sq, for k=1 and k=K_B
lambda_k.sq.full =c(lambda0.sq, lambda_k.sq, lambda0.sq)

############################
## Matrices to hold draws ## 
############################
MCMC.BASOFR = structure(list(B.post=array(NA, c(S, dim(X)[2])),
                             sigma.2.post = array(NA, c(S,1)),
                             lambda.sq.post = array(NA, c(S, K)), 
                             ##  DHS paramaters
                             mu.post = rep(NA, S), # unconditional mean of the AR 1 process
                             phi.post = rep(NA,S), # persistence parameter of the AR 1 process
                             h_sd.post = array(NA, c(S,(K-2))))) # innovation terms of the AR 1 process

#####################
## Begin MCMC Loop ##
#####################
D.matrix = rbind(c(1,rep(0,ncol(D_2)-1)), D_2, c(rep(0,ncol(D_2)-1),1))

t0 = Sys.time()
for (s in 1:S) {
  ## update sigma.2
  a.new = a0 + n/2
  b.new = b0 + sum((y- crossprod(t(X),B))^2)/2
  sigma.2=1/rgamma(1, shape=a.new, rate=b.new)
  
  ## update B
  Q.B = crossprod(X)/sigma.2 +rbind(0,cbind(0,crossprod(D.matrix, 1/lambda_k.sq.full *D.matrix))) 
  l.B = crossprod(X,y)/sigma.2
  ch.QB = chol(Q.B)
  B = backsolve(ch.QB, forwardsolve(t(ch.QB), l.B) + rnorm(length(l.B)))
  
  ## update {lambda_k.sq} and all the DHS parameters (mu, phi and h_sd)
  omega = as.matrix(tcrossprod(D_2, t(as.matrix(B[-1]))))
  sample = sampleDSP(omega = omega,evolParams = evolParams)
  h_k  = sample$ht # log (lambda_k.sq)
  lambda_k.sq = as.numeric(exp(h_k))
  mu = sample$dhs_mean; phi = sample$dhs_phi
  sigma_eta_t = sample$sigma_eta_t; sigma_eta_0 = sample$sigma_eta_0
  
  evolParams = list(ht = h_k,
                    dhs_mean = mu,
                    dhs_phi = phi,
                    sigma_eta_t = sigma_eta_t,
                    sigma_eta_0 = sigma_eta_0)
  
  ## sample lambda_0.sq
  a.lambda.new = a.lambda +1
  b.lambda.new = b.lambda + (B[2]^2 + B[length(B)]^2)/2
  lambda0.sq=1/rgamma(1, a.lambda.new, b.lambda.new)
  
  ## update the lambda_k.sq.full
  lambda_k.sq.full = c(lambda0.sq, lambda_k.sq, lambda0.sq)
  
  MCMC.BASOFR$B.post[s,] = B
  MCMC.BASOFR$sigma.2.post[s] = sigma.2
  MCMC.BASOFR$lambda.sq.post[s,]= lambda_k.sq.full
  MCMC.BASOFR$mu.post[s] = mu;
  MCMC.BASOFR$phi.post[s] = phi;
  MCMC.BASOFR$h_sd.post[s,] = c(sigma_eta_0, sigma_eta_t)
}
print(paste0("MCMC running time is :",Sys.time() - t0))

## To replicate the result from Figure 3,
## we saved the results for n=500, SNR=5, smoothness = "smooth" as:
# save(MCMC.BASOFR, file = "MCMC_BASOFR_n500_SNR5_smooth.RDATA")

#######################
## Visualize the fit ##
#######################

## discard the burn-in, compute the posterior distribution of the regression function estimate
beta.post.sim = tcrossprod(tail(MCMC.BASOFR$B.post,iter)[,-1],basismat)

## posterior mean estimate of the regression function
beta.pm  = apply(beta.post.sim, 2, mean)

## 95% and 50% credible bands estimate
i=95
beta.low95 = apply(beta.post.sim, 2, function(x) quantile(x, (1-i/100)/2))
beta.high95 = apply(beta.post.sim, 2, function(x) quantile(x, 1-(1-i/100)/2))
i=50
beta.low50 = apply(beta.post.sim, 2, function(x) quantile(x, (1-i/100)/2))
beta.high50 = apply(beta.post.sim, 2, function(x) quantile(x, 1-(1-i/100)/2))

## Create the plot
title = paste0("BASOFR, n=",n, ", SNR=", STN, ", ",smoothness )
ggplot() + 
  coord_cartesian(ylim = c(-3.25, 2)) + 
  geom_line(aes(x=x,y=beta), col="red", size=1.5)+
  geom_line(aes(x=x, y=beta.pm), size=1.5) + 
  geom_ribbon(aes(x=x, y=NULL, ymin = beta.low95, ymax=beta.high95), fill="blue", alpha=I(1/5))+
  geom_ribbon(aes(x=x, y=NULL, ymin = beta.low50, ymax= beta.high50), fill="blue", alpha=I(2/5))+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size= 20, face= "bold")) + 
  ylab(TeX("$\\beta(t)$")) + xlab(TeX("$\\t"))+
  ggtitle(title)
