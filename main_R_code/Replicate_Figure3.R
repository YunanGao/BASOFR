################################################
## Code to Replicate the result from Figrue 3 ##
################################################

###################
## Preliminaries ##
###################

##Load Libraries
library(ggplot2)
library(fda)
library(MASS)
library(latex2exp)

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
set.seed(1)
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

##################################################
## MCMC settings for the local P-spline  method ##
##################################################
burn <- 10000
iter <- 10000

## total iterations
S = burn + iter 

######################
## Prior parameters ##
######################

## Gamma prior on the SOFR (1/{\sigma^2}) term
a0 = 0.01; b0=0.01

## Gamma prior on {\lambda_k}^{-2}, for k=2,...,(K_B-1)
alpha.k=rep(0.01,K-2); beta.k=rep(0.01,K-2)


#####################
## Starting values ##
#####################

## the intercept term and the B-spline coefficients of the regression function
B=crossprod(t(solve(crossprod(X) + diag(10^(-5), nrow=(K+1)))),crossprod(X, y))

## lambda_k.sq, for k=2,...,(K_B-1)
lambda_k.sq = rep(1, K-2)

############################
## Matrices to hold draws ## 
############################
MCMC.LPSOFR = structure(list(B.post=array(NA, c(S, dim(X)[2])),
                             sigma.2.post = array(NA, c(S,1)),
                             lambda.sq.post = array(NA, c(S, K-2)))) 

#####################
## Begin MCMC Loop ##
#####################

D_2 = diff(diag(K), diff=2)

for (s in 1:S) {
  # update sigma.2
  a.new = a0 + n/2
  b.new = b0 + sum((y- crossprod(t(X),B))^2)/2
  sigma.2=1/rgamma(1, shape=a.new, rate=b.new)
  
  # update B
  Q.B = crossprod(X)/sigma.2 + rbind(0,cbind(0,crossprod(D_2, 1/lambda_k.sq * D_2))) 
  l.B = crossprod(X,y)/sigma.2
  ch.QB = chol(Q.B)
  B = backsolve(ch.QB, forwardsolve(t(ch.QB), l.B) + rnorm(length(l.B)))
  
  # update {lambda_k.sq}, k=2,...,(K_B-1)
  omega = as.matrix(tcrossprod(D_2, t(as.matrix(B[-1]))))
  alpha.k.new = alpha.k  + 1/2
  beta.k.new = beta.k + 0.5*omega^2
  lambda_k.sq = 1/rgamma(K-2, shape = alpha.k.new, rate = beta.k.new)
  
  MCMC.LPSOFR$B.post[s,] = B
  MCMC.LPSOFR$sigma.2.post[s] = sigma.2
  MCMC.LPSOFR$lambda.sq.post[s,]= lambda_k.sq
}


########################
## Reproduce Figure 3 ##
########################

## Load the BASOFR MCMC result (This is produced by the fit_BASOFR_Sim.R code)
load("MCMC_BASOFR_n500_SNR5_smooth.RDATA")

## discard the burn-in, compute the posterior distribution of the regression function estimate 
beta.post.sim.BASOFR = tcrossprod(tail(MCMC.BASOFR$B.post,iter)[,-1],basismat)
beta.post.sim.LPSOFR = tcrossprod(tail(MCMC.LPSOFR$B.post,iter)[,-1],basismat)

## posterior mean estimate of the regression function 
beta.pm.BASOFR  = apply(beta.post.sim.BASOFR, 2, mean)
beta.pm.LPSOFR  = apply(beta.post.sim.LPSOFR, 2, mean)

## 95% and 50% credible bands estimate 
i=95
beta.low95.BASOFR = apply(beta.post.sim.BASOFR, 2, function(x) quantile(x, (1-i/100)/2))
beta.high95.BASOFR = apply(beta.post.sim.BASOFR, 2, function(x) quantile(x, 1-(1-i/100)/2))
beta.low95.LPSOFR = apply(beta.post.sim.LPSOFR, 2, function(x) quantile(x, (1-i/100)/2))
beta.high95.LPSOFR = apply(beta.post.sim.LPSOFR, 2, function(x) quantile(x, 1-(1-i/100)/2))
i=50
beta.low50.BASOFR = apply(beta.post.sim.BASOFR, 2, function(x) quantile(x, (1-i/100)/2))
beta.high50.BASOFR = apply(beta.post.sim.BASOFR, 2, function(x) quantile(x, 1-(1-i/100)/2))
beta.low50.LPSOFR = apply(beta.post.sim.LPSOFR, 2, function(x) quantile(x, (1-i/100)/2))
beta.high50.LPSOFR = apply(beta.post.sim.LPSOFR, 2, function(x) quantile(x, 1-(1-i/100)/2))

## Create the fitting plot for the BASOFR method
title.BASOFR = paste0("BASOFR, n=",n, ", SNR=", STN, ", ",smoothness )
ggplot() + 
  coord_cartesian(ylim = c(-3.25, 2)) + 
  geom_line(aes(x=x,y=beta), col="red", size=1.5)+
  geom_line(aes(x=x, y=beta.pm.BASOFR), size=1.5) + 
  geom_ribbon(aes(x=x, y=NULL, ymin = beta.low95.BASOFR, ymax=beta.high95.BASOFR), fill="blue", alpha=I(1/5))+
  geom_ribbon(aes(x=x, y=NULL, ymin = beta.low50.BASOFR, ymax= beta.high50.BASOFR), fill="blue", alpha=I(2/5))+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size= 20, face= "bold")) + 
  ylab(TeX("$\\beta(t)$")) + xlab(TeX("$\\t"))+
  ggtitle(title.BASOFR)

## Create the fitting plot for the local P-spline method
title.LPSOFR = paste0("Local P-spline, n=",n, ", SNR=", STN, ", ",smoothness )
ggplot() + 
  coord_cartesian(ylim = c(-4, 4)) + 
  geom_line(aes(x=x,y=beta), col="red", size=1.5)+
  geom_line(aes(x=x, y=beta.pm.LPSOFR), size=1.5) + 
  geom_ribbon(aes(x=x, y=NULL, ymin = beta.low95.LPSOFR, ymax=beta.high95.LPSOFR), fill="blue", alpha=I(1/5))+
  geom_ribbon(aes(x=x, y=NULL, ymin = beta.low50.LPSOFR, ymax= beta.high50.LPSOFR), fill="blue", alpha=I(2/5))+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size= 20, face= "bold")) + 
  ylab(TeX("$\\beta(t)$")) + xlab(TeX("$\\t"))+
  ggtitle(title.LPSOFR)
