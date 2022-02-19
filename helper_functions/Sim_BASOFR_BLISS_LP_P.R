Sim_BASOFR_LP_P_BLISS = function(N=50, #N is the simulation number
                         n, STN, smoothness){
  
  ## Are the functional predictors smooth or nonsmooth?
  smoothness = "smooth"
  if(smoothness == "smooth"){l = 0.1}
  if(smoothness == "nonsmooth"){l=0.01}
  
  ## Set true regression function \beta(t)
  beta = 8 * (2+exp(20-60*x) + exp(60*x-20))^(-1) - 12*(2 + exp(40-60*x) + exp(60*x-40))^(-1)
  
  ## Matrix to store the simulation result
  if(n > 100){ # if n is larger than 100, the competitors do not include BLISS, since BLISS would take too long time
    simulation.result = matrix(NA, nrow=0, ncol=9)
    colnames(simulation.result) = c("L2_BASOFR", "L2_LP", "L2_P",
                                    "cov_BASOFR", "cov_LP", "cov_P",
                                    "area_BASOFR", "area_LP", "area_P")
  }else{
    simulation.result = matrix(NA, nrow=0, ncol=12)
    colnames(simulation.result) = c("L2_BASOFR", "L2_LP", "L2_P", "L2_BLISS",
                                    "cov_BASOFR", "cov_LP", "cov_P", "cov_BLISS",
                                    "area_BASOFR", "area_LP", "area_P", "area_BLISS")
  }
  
  
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
  
  ## estimate the J matrix
  J = eval.penalty(basisobj= basis, Lfdobj=0, rng=c(0,1))
  
  ## Start looping over from simulation 1 to simulation N
  for (N_i in 1:N) {
    set.seed(N_i)
    beta = 8 * (2+exp(20-60*x) + exp(60*x-20))^(-1) - 12*(2 + exp(40-60*x) + exp(60*x-40))^(-1)
    
    ## Simulate the functional observations from a Gaussian process
    x_diff = abs(matrix(rep(x,obs.x), nrow=obs.x)-t(matrix(rep(x,obs.x),nrow=obs.x)))
    sigma_x = 0.7
    Kernel = sigma_x^2 * exp(-x_diff^2/(2*l^2))
    mu = rep(0,obs.x)
    x.sim = mvrnorm(n=n, mu=mu, Sigma=Kernel) 
    
    ## Simulate the scalar responses 
    noise.sd = sd(tcrossprod(x.sim, t(beta))/(obs.x-1))/STN
    y = tcrossprod(x.sim, t(beta))/(obs.x-1) + rnorm(n, mean=0,sd=noise.sd)
    
    ## estimate the OLS b-spline basis coefficients
    X.star = t(tcrossprod(tcrossprod(solve(crossprod(basismat)), basismat), x.sim))
    
    ## the predictors for the equivalent multiple linear regression
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
    
    #########################
    ## MCMC for the BASOFR ##
    #########################
    
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
    
    ########################################
    ## MCMC for the Local P-spline method ##
    ########################################
    
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
    
    ########################################
    ## MCMC for the P-spline method ##
    ########################################
    
    ######################
    ## Prior parameters ##
    ######################
    
    ## Gamma prior on the SOFR (1/{\sigma^2}) term
    a0 = 0.01; b0=0.01
    
    ## Gamma prior on {\lambda}^{-2}, where \lambda^2 is the global prior variance on the second-differenced B-spline coefficients
    alpha.0 = 0.01; beta.0 = 0.01
    
    #####################
    ## Starting values ##
    #####################
    ## the intercept term and the B-spline coefficients of the regression function
    B=crossprod(t(solve(crossprod(X) + diag(10^(-5), nrow=(K+1)))),crossprod(X, y))
    
    ## global scale parameter
    lambda_k.sq = rep(1, K-2)
    
    ############################
    ## Matrices to hold draws ## 
    ############################
    MCMC.PSOFR = structure(list(B.post=array(NA, c(S, dim(X)[2])),
                             sigma.2.post = array(NA, c(S,1)),
                             lambda.sq.post = array(NA, c(S, (K-2)))))
    
    
    #####################
    ## Begin MCMC Loop ##
    #####################
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
      
      # update {lambda_0.sq}
      omega = as.matrix(tcrossprod(D_2, t(as.matrix(B[-1]))))
      alpha.new = alpha.0 + length(omega)/2
      beta.new = sum(omega^2)/2 + beta.0
      lambda.sq = 1/rgamma(1, shape=alpha.new, rate = beta.new)
      lambda_k.sq = rep(lambda.sq, length(omega))
      
      MCMC.PSOFR$B.post[s,] = B
      MCMC.PSOFR$sigma.2.post[s] = sigma.2
      MCMC.PSOFR$lambda.sq.post[s,]= lambda_k.sq
    }
    
    ## discard the burn-in, compute the posterior distribution of the regression function estimate 
    beta.post.sim.BASOFR = tcrossprod(tail(MCMC.BASOFR$B.post,iter)[,-1],basismat)
    beta.post.sim.LPSOFR = tcrossprod(tail(MCMC.LPSOFR$B.post,iter)[,-1],basismat)
    beta.post.sim.PSOFR = tcrossprod(tail(MCMC.PSOFR$B.post,iter)[,-1],basismat)
    
    ## posterior mean estimate of the regression function 
    beta.pm.BASOFR  = apply(beta.post.sim.BASOFR, 2, mean)
    beta.pm.LPSOFR  = apply(beta.post.sim.LPSOFR, 2, mean)
    beta.pm.PSOFR = apply(beta.post.sim.PSOFR, 2, mean)
    
    ## 95% and 50% credible bands estimate 
    i=95
    beta.low95.BASOFR = apply(beta.post.sim.BASOFR, 2, function(x) quantile(x, (1-i/100)/2))
    beta.high95.BASOFR = apply(beta.post.sim.BASOFR, 2, function(x) quantile(x, 1-(1-i/100)/2))
    beta.low95.LPSOFR = apply(beta.post.sim.LPSOFR, 2, function(x) quantile(x, (1-i/100)/2))
    beta.high95.LPSOFR = apply(beta.post.sim.LPSOFR, 2, function(x) quantile(x, 1-(1-i/100)/2))
    beta.low95.PSOFR = apply(beta.post.sim.PSOFR, 2, function(x) quantile(x, (1-i/100)/2))
    beta.high95.PSOFR = apply(beta.post.sim.PSOFR, 2, function(x) quantile(x, 1-(1-i/100)/2))
    i=50
    beta.low50.BASOFR = apply(beta.post.sim.BASOFR, 2, function(x) quantile(x, (1-i/100)/2))
    beta.high50.BASOFR = apply(beta.post.sim.BASOFR, 2, function(x) quantile(x, 1-(1-i/100)/2))
    beta.low50.LPSOFR = apply(beta.post.sim.LPSOFR, 2, function(x) quantile(x, (1-i/100)/2))
    beta.high50.LPSOFR = apply(beta.post.sim.LPSOFR, 2, function(x) quantile(x, 1-(1-i/100)/2))
    beta.low50.PSOFR = apply(beta.post.sim.PSOFR, 2, function(x) quantile(x, (1-i/100)/2))
    beta.high50.PSOFR = apply(beta.post.sim.PSOFR, 2, function(x) quantile(x, 1-(1-i/100)/2))
    
    ## L2 error
    l2.error.BASOFR= sum((beta - beta.pm.BASOFR)^2/(obs.x-1))
    l2.error.LPSOFR= sum((beta - beta.pm.LPSOFR)^2/(obs.x-1))
    l2.error.PSOFR= sum((beta - beta.pm.PSOFR)^2/(obs.x-1))
    
    ## 95% CI area
    CI95.area.BASOFR = sum((beta.high95.BASOFR -beta.low50.BASOFR)/(obs.x-1))
    CI95.area.LPSOFR = sum((beta.high95.LPSOFR -beta.low50.LPSOFR)/(obs.x-1))
    CI95.area.PSOFR = sum((beta.high95.PSOFR -beta.low50.PSOFR)/(obs.x-1))
    
    ## Coverage of the 95% CI 
    cov.BASOFR = mean(beta.low95.BASOFR <= beta & beta.high95.BASOFR >= beta)
    cov.LPSOFR = mean(beta.low95.LPSOFR <= beta & beta.high95.LPSOFR >= beta)
    cov.PSOFR = mean(beta.low95.PSOFR <= beta & beta.high95.PSOFR >= beta)
    
    result.N_i = c(l2.error.BASOFR, l2.error.LPSOFR, l2.error.PSOFR, 
                   cov.BASOFR, cov.LPSOFR, cov.PSOFR, 
                   CI95.area.BASOFR, CI95.area.LPSOFR, CI95.area.PSOFR)
    
    if(n <= 100){
      ## Run Bliss
      
      ## Create the data in the form fit_Bliss require
      data = list()
      x.grid = seq(0,1,length.out=obs.x)
      data$Q = 1; data$y = as.numeric(y); data$x = list(x.sim); data$grids = list(x.grid)
      
      ## BIC model selection
      iter = 1e3
      Ks = 1:5
      BIC.bliss = BIC_model_choice(Ks, iter, data)
      K.value = which(BIC.bliss == min(BIC.bliss))
      
      ## Param for the Gibbs
      param <- list(
        iter = 1e3,
        burnin = 5e2,
        K=K.value
      )
      
      ## fit the BLISS model
      fit_Bliss <- fit_Bliss(data=data, param=param, verbose=TRUE)
      
      ## the posterior distribution of the regression function estimate
      beta.post.bliss = fit_Bliss$beta_posterior_density[[1]]$new_beta_sample
      
      ## 95% and 50% credible bands estimate
      beta.bliss.high95 = apply(beta.post.bliss, 2, function(x) quantile(x,0.975))
      beta.bliss.low95 = apply(beta.post.bliss, 2, function(x) quantile(x,0.025))
      beta.bliss.high50 = apply(beta.post.bliss, 2, function(x) quantile(x,0.75))
      beta.bliss.high50 = apply(beta.post.bliss, 2, function(x) quantile(x,0.25))
      
      ## L2 error of Bliss method
      l2.error.BLISS = sum((beta-fit_Bliss$Smooth_estimate[[1]])^2/(obs.x-1))
      
      ## 95% CI area of the Bliss method
      CI95.area.BLISS = sum(beta.bliss.high95 - beta.bliss.low95)/(obs.x-1)
      
      ## Coverage of the Bliss method      
      beta[which(abs(beta)<=5e-4)]=0
      cov.BLISS = mean(beta.bliss.low95 <= beta & beta.bliss.high95 >= beta)

      result.N_i = c(l2.error.BASOFR, l2.error.LPSOFR, l2.error.PSOFR, l2.error.BLISS,
                     cov.BASOFR, cov.LPSOFR, cov.PSOFR, cov.BLISS,
                     CI95.area.BASOFR, CI95.area.LPSOFR, CI95.area.PSOFR, CI95.area.BLISS)
    }
    
    simulation.result = rbind(simulation.result, result.N_i)
    
  }
  print("Job is done")
  return(simulation.result)
  
}



