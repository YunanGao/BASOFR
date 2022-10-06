###########################################
## Decision analysis tools for selecting ##
##      the critical windows             ##
###########################################

# We use this code to show an example of using the proposed decision analysis strategy to
# identify critical windows of susceptibility from a Bayesian SOFR. This code contains the following contents:

# 1. estimate an BASOFR model
# 2. extract the critical windows from the model using the proposed decision analysis tools.
# 3. evalute the critical window selection using true positive rate (TPR) and true negative rate (TNR)
# and compare the results between the decision analysis approach and the credible interval-based approach. 


###################
## Preliminaries ##
###################

##Load Libraries
library(MASS)
library(fda)
library(ggplot2)
library(coda)
library(stochvol)
library(BayesLogit)
library(genlasso)
library(latex2exp)
library(RColorBrewer)

## Source in Helper Functions
## Users need to set the R working directory as the path of the main_R_code file
setwd('/Users/yunangao/Downloads/BASOFR-main 4/main_R_code')
source("../helper_functions/dhs_sampling.R")

## Set sample size, signal-to-noise ratio (SNR)
n = 100000; STN = 0.05

## Set true regression function \beta(t)
beta = -0.02* as.numeric ( x <=0.333) + 0.02*as.numeric(x> 0.333 & x<=0.67) - 0.02*as.numeric(x>=0.67 & x<=1)

#######################
## simulate the data ##
#######################

# This R function below generate functional covariates with seasonality patterns
FunctionalPredictors= function(n, l, x=seq(0,1,length.out = 101), sigma_x = 0.7, seasonality=TRUE,  P=365/(295-1), A=1 ){
  ## phase parameters
  if(seasonality == TRUE){
    phi.phase = as.matrix(runif(n, min=0, max=2*pi))
    
    ## mean functions for all n subjects
    mu.x.matrix = apply(phi.phase, 1, function(phi) A*sin((2*pi*x/P) + phi))
  }
  
  
  ## Kernels
  obs.x = length(x)
  x_diff= abs(matrix(rep(x, obs.x), nrow=obs.x) - t(matrix(rep(x, obs.x), nrow=obs.x)))
  Kernel = sigma_x^2 * exp(-x_diff^2/(2*l^2))
  
  
  ## functional predictors
  if(seasonality == TRUE){
    x.sim = t(apply(mu.x.matrix, 2, function(mu.x) mvrnorm(n=1, mu=mu.x, Sigma=Kernel)))
    return(x.sim)
  }else{
    mu = rep(0, obs.x)
    x.sim = mvrnorm(n=n, mu=mu, Sigma=Kernel)
    return(x.sim)}
  
}


## this code simulate functional covariates
set.seed(1) 
x.sim = FunctionalPredictors(n=100000, l=0.001, seasonality = TRUE,  P=365/((295-1)*2))

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
burn <- 1000
iter <- 1000

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

#######################
## Visualize the fit ##
#######################

beta.post.sim.BASOFR = tcrossprod(tail(MCMC.BASOFR$B.post,iter)[,-1],basismat)
beta.pm.BASOFR  = apply(beta.post.sim.BASOFR, 2, mean)
i=95
beta.low95.BASOFR = apply(beta.post.sim.BASOFR, 2, function(x) quantile(x, (1-i/100)/2))
beta.high95.BASOFR = apply(beta.post.sim.BASOFR, 2, function(x) quantile(x, 1-(1-i/100)/2))
i=50
beta.low50.BASOFR = apply(beta.post.sim.BASOFR, 2, function(x) quantile(x, (1-i/100)/2))
beta.high50.BASOFR = apply(beta.post.sim.BASOFR, 2, function(x) quantile(x, 1-(1-i/100)/2))

BASOFR.plot = ggplot() + 
  coord_cartesian(ylim = c(-0.08, 0.06)) + 
  geom_line(aes(x=x, y=beta), col="red", size=1.5) + 
  geom_line(aes(x=x, y=beta.pm.BASOFR),  size=2) + 
  geom_ribbon(aes(x=x, y=NULL, ymin = beta.low95.BASOFR, ymax = beta.high95.BASOFR), fill="blue", alpha = I(1/5))+
  geom_ribbon(aes(x=x, y=NULL, ymin = beta.low50.BASOFR, ymax = beta.high50.BASOFR), fill="blue", alpha = I(2/5))+
  theme(axis.text = element_text(size=20),
          axis.title = element_text(size=20 ,face = "bold")) + 
  ylab(TeX("$\\beta(t)")) + xlab(TeX("$\\t"))

BASOFR.plot 

# the BASOFR posterior mean and credible interval offer limited ability to descrbe the true regression function

#############################
## Decision Analysis tools ##
#############################

## Suppose delta(t) is locally constant over a partition, and we calculate the aggregated trajectories of the X(t) over these subdomains.
partition.K = seq(0,1,0.05)
K = length(partition.K) -1
x.subdomain.aver.K = t(apply(x.sim, 1, function(x)
  rowSums(matrix(x[-1], ncol=(1/K)/(1/(obs.x-1)), byrow = TRUE))*0.01 + 
    c(x[1]*0.01,rep(0,(K-1)))))
intercept_x.subdomain.aver.K = cbind(1, x.subdomain.aver.K)

## Obtain the predictive posterior distribution of {y_i-\mu}
B.post = MCMC.BASOFR$B.post
sigma.2.post = MCMC.BASOFR$sigma.2.post
beta.post.sim = tail(beta.post.sim.BASOFR,iter)
beta.pm = beta.pm.BASOFR
intercept.post = tail(B.post[,1], iter)
intercept.pm = mean(intercept.post)
sd.post = sqrt(tail(sigma.2.post,iter))
y.predictive.posterior = tcrossprod(x.sim, beta.post.sim)/(obs.x-1)+ 
  apply(sd.post, 1, function(x) rnorm(n, mean=0, sd=x))

## the posterior mean of {y_i - \mu} 
h.bar = x.sim%*%beta.pm/(obs.x-1)

## use the genlasso R package to obtain multiple locally constant estimates of different complexity (i.e. solve eq (10)) in the manuscript)  
d=1
D = rbind(0,cbind(0, diff(diag(K), diff=d))) 
out = fusedlasso(y=h.bar, X=intercept_x.subdomain.aver.K, D=D)

## Collect the changing points of the locally constant approximations
num.cp = apply(diff(diag(K), diff=1) %*% out$beta[-1,], 2, function(x)
  sum(x!=0))
cp.location=(diff(diag(K), diff=1) %*% out$beta[-1,] != 0) * 1
model.sum = unique(as.data.frame(cbind(num.cp, t(cp.location))))
model.sum = model.sum[order(model.sum$num.cp),]

## Matrices to hold empirical MSE and predictive MSE
MSE.predictive.insample = array(NA, dim=c(dim(model.sum)[1],iter))
MSE.insample = rep(NA, dim(model.sum)[1])

## Obtain the predictive and empirical MSE for each locally constant estimates and 
par(mfrow=c(2,2))
for (m in 1:dim(model.sum)[1]) {
  model = model.sum[m,]
  
  ## obtain an unbiased fit for each locally constant approximation
  domain.m = matrix(c(0, 
                      sort(c(which(model[-1]!=0)/K, which(model[-1]!=0)/K+0.01 )),
                      1) *(obs.x-1) + 1,
                    ncol=2, byrow=T)
  x.m = apply(domain.m, 1, function(domain)
    apply(x.sim, 1, function(x) sum(x[domain[1]:domain[2]]))/obs.x) 
  X.m = cbind(1,x.m)
  lm.m = lm(h.bar~x.m)
  times = round(domain.m[,2]-domain.m[,1]+1)
  beta.m = rep(c(lm.m$coefficients[-1]), times)
  if(length(beta.m) != obs.x){
    print("warning! approx beta length != obs x length")
  }
  
  ## Visualize the locally constant estimate
  main = paste0("model", m, " with npc=", model[1])
  plot(x, beta.m, type="l", lwd=4,xlab="t", ylab="beta(t)",main=main, ylim = c(min(beta.m, beta.pm, beta),
                                                                               max(beta.m, beta.pm, beta)))
  lines(x,beta, type="b",col="red")
  lines(x, beta.pm, col="blue",lty=6,lwd=3)
  legend("bottomright", legend=c("BAM", "BASOFR", "truth"), 
         col=c("black", "blue", "red"), lty=c(6,6,6), cex=0.8, lwd=c(4,3,3),
         bg='transparent')
  
  ## predictive MSE and empirical MSE
  MSE.predictive.insample[m, ]= apply(y.predictive.posterior, 2, function(x) 
    mean((x-lm.m$fitted.values)^2)) 
  MSE.insample[m] = mean((y-intercept.pm - lm.m$fitted.values)^2)
}

## Epsilon level: sparsest model w/prob. eps_level of matching predictions for best model 
eps_level = 0.1

## Best model by in-sample MSE:
ell_min_in = which.min(MSE.insample)

## Percent difference in predictive MSE relative to the "best" model:
post.dmse = apply(t(MSE.predictive.insample), 2 , function(pmse)
  100*(pmse- t(MSE.predictive.insample)[,ell_min_in])/t(MSE.predictive.insample)[,ell_min_in])

## Percent difference in MSE relative to the "best" model:
dmse = 100*(MSE.insample - MSE.insample[ell_min_in])/MSE.insample[ell_min_in]

## 100(1 - 2*eps_level)% 
ci_dmse = t(apply(post.dmse, 2, quantile, c(eps_level, 1 - eps_level)))

## one of the Smallest model(s) for which the lower CI includes zero:
ell_eps_in = min(which(colMeans(post.dmse <= 0) >= eps_level))

## num of Cps in the ell_eps_in-th model
ncp.ell_eps_in = model.sum$num.cp[ell_eps_in]

## all the models' index with the same num of Cps
same.size.index.in = which(model.sum$num.cp == ncp.ell_eps_in) 

## the set of models we are interested in
models.int.in = same.size.index.in[same.size.index.in >= ncp.ell_eps_in]

## the smallest model(s) for which the lower CI includes zero:
ell_eps_in_set = intersect(which(colMeans(post.dmse <= 0) >= eps_level),
                           models.int.in)
par(mfrow=c(1,1))
model.index = 1:dim(model.sum)[1]

###################################
# Selecting the acceptable family #                  
################################### 
par(mai=c(1,1.5,1,0.5))
plot(model.index,
     colMeans(post.dmse),type='p', ylim = range(ci_dmse, 0), lwd=5,
     xlab = 'Model Index', ylab = 'Difference in predictive loss(%)', main = '',
     cex.lab = 2, cex.axis=2)
abline(h = 0, lwd=3, lty=6)
abline(v = ell_eps_in_set, lwd = 3, col="darkgrey")
abline(v = ell_min_in, lwd=3, col="lightgray", lty=3)
arrows(model.index[-ell_min_in], ci_dmse[-ell_min_in,1],
       model.index[-ell_min_in], ci_dmse[-ell_min_in,2],
       length=0.05, angle=90, code=3, lwd=3)
lines(model.index, dmse, type='p', lwd=7,col='gray', pch=4, cex=3)

# All the model(s) on the right hand side of the dark grey model belongs to the acceptable family (epsilon),
# and they all have competitive predictive performanaces. The dark grey model is the simplest member among the acceptable family

#############################################################
## Visualize the simplest member in the acceptable family  ##
#############################################################

m=ell_eps_in_set
model = model.sum[m,]
domain.m = matrix(c(0, 
                    sort(c(which(model[-1]!=0)/K, which(model[-1]!=0)/K+0.01 )),
                    1) *(obs.x-1) + 1,
                  ncol=2, byrow=T)

# Define the new covariate (X.m) of model m, dim(X.m) = n x (ncps.m + 2)
x.m = apply(domain.m, 1, function(domain)
  apply(x.sim, 1, function(x) sum(x[domain[1]:domain[2]]))/obs.x) 
X.m = cbind(1,x.m)

lm.m = lm(h.bar~x.m)
times = round(domain.m[,2]-domain.m[,1]+1)
beta.DA = rep(c(lm.m$coefficients[-1]), times)

ggplot() + 
  geom_hline(yintercept=0, col= 'white', size=2)+
  coord_cartesian(ylim = c(-0.08, 0.06)) + 
  geom_line(aes(x=x, y=beta), col="red", size=1.5) + 
  geom_line(aes(x=x, y=beta.pm.BASOFR),  size=2) + 
  geom_line(aes(x=x, y=beta.DA), size=2, col='darkgreen')+
  geom_ribbon(aes(x=x, y=NULL, ymin = beta.low95.BASOFR, ymax = beta.high95.BASOFR), fill="blue", alpha = I(1/5))+
  geom_ribbon(aes(x=x, y=NULL, ymin = beta.low50.BASOFR, ymax = beta.high50.BASOFR), fill="blue", alpha = I(2/5))

# Though the posterior mean and credible intervals offer limited ability to describe the true regression function,
# the proposed decision analysis approach adequately recovers the truth. 
        
#############################################################
## Compare the critical window selection using the proposed decision  ##
## analysis approch v.s. credible interval-based selection
#############################################################


## probability of positive / negative critical windows in the true data
pr.positive = sum(beta>0)/length(beta)
pr.negative = 1-pr.positive

## True positive rate (TPR)
true.positive.CI = sum(beta.low95.BASOFR>0 & beta > 0)/length(beta)
true.positive.DA = sum(beta.DA >0 & beta>0)/length(beta)
TPR.CI = true.positive.CI/pr.positive
TPR.DA = true.positive.DA/pr.positive

## True negative rate (TNR)
true.negative.CI=  sum(beta.high95.BASOFR < 0 & beta<0)/length(beta)
true.negative.DA = sum(beta.DA <0 & beta <0)/length(beta)
TNR.CI = true.negative.CI/pr.negative
TNR.DA = true.negative.DA/pr.negative

## L2-error (posterior mean v.s. locally constant estimate)

# l2-error of the posterior mean
(l2.error.BASOFR = sum((beta - beta.pm.BASOFR)^2/(obs.x-1)))
# l2-error of the locally constant estimate
(l2.error.DA = sum((beta - beta.DA)^2/(obs.x-1)))

df.result = data.frame(matrix(NA, nrow=2, ncol=2))
colnames(df.result) = c('TPR', 'TNR')
rownames(df.result) = c('CI', 'DA')

df.result[1,] = c(TPR.CI, TNR.CI)
df.result[2,] = c(TPR.DA, TNR.DA)
df.result




