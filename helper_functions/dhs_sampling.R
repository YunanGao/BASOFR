# These functions sample the DHS compoenents.
# These functions come from the R code of the Kowal, Daniel R., David S. Matteson, and David Ruppert.
# "Dynamic shrinkage processes." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 81, no. 4 (2019): 781-804.
# Link: https://github.com/drkowal/dsp/blob/master/R/helper_functions.R
#----------------------------------------------------------------------------
#' Initialize the evolution error variance parameters
#'
#' Compute initial values for evolution error variance parameters under the dynamic horseshoe prior
#'
#' @param omega \code{T x p} matrix of evolution errors
#' @return List of relevant components: the \code{T x p} evolution error SD \code{sigma_wt},
#' the \code{T x p} log-volatility \code{ht}, the \code{p x 1} log-vol unconditional mean(s) \code{dhs_mean},
#' the \code{p x 1} log-vol AR(1) coefficient(s) \code{dhs_phi},
#' the \code{T x p} log-vol innovation SD \code{sigma_eta_t} from the PG priors,
#' the \code{p x 1} initial log-vol SD \code{sigma_eta_0},
#' and the mean of log-vol means \code{dhs_mean0} (relevant when \code{p > 1})
#' @export
#' 
initDHS = function(omega){
  
  # "Local" number of time points
  omega = as.matrix(omega)
  n = nrow(omega); p = ncol(omega)
  
  # Initialize the log-volatilities:
  ht = log(omega^2 + 0.0001)
  
  # Initialize the AR(1) model to obtain unconditional mean and AR(1) coefficient
  arCoefs = apply(ht, 2, function(x){
    params = try(arima(x, c(1,0,0))$coef, silent = TRUE); if(class(params) == "try-error") params = c(0.8, mean(x)/(1 - 0.8))
    params
  })
  dhs_mean = arCoefs[2,]; dhs_phi = arCoefs[1,]; dhs_mean0 = mean(dhs_mean)
  
  # Initialize the SD of log-vol innovations simply using the expectation:
  sigma_eta_t = matrix(pi, nr = n-1, nc = p)
  sigma_eta_0 = rep(pi, p) # Initial value
  
  # Evolution error SD:
  sigma_wt = exp(ht/2)
  
  list(sigma_wt = sigma_wt, ht = ht, dhs_mean = dhs_mean, dhs_phi = dhs_phi, sigma_eta_t = sigma_eta_t, sigma_eta_0 = sigma_eta_0, dhs_mean0 = dhs_mean0)
}

#' Sample components from a discrete mixture of normals
#'
#' Sample Z from 1,2,...,k, with P(Z=i) proportional to q[i]N(mu[i],sig2[i]).
#'
#' @param y vector of data
#' @param mu vector of component means
#' @param sig vector of component standard deviations
#' @param q vector of component weights
#' @return Sample from {1,...,k}
#----------------------------------------------------------------------------
ncind = function(y,mu,sig,q){
  sample(1:length(q),
         size = 1,
         prob = q*dnorm(y,mu,sig))
}
#----------------------------------------------------------------------------
#' Sample the AR(1) coefficient(s)
#'
#' Compute one draw of the AR(1) coefficient in a model with Gaussian innovations
#' and time-dependent innovation variances. In particular, we use the sampler for the
#' log-volatility AR(1) process with the parameter-expanded Polya-Gamma sampler. The sampler also applies
#' to a multivariate case with independent components.
#'
#' @param h_yc the \code{T x p} matrix of centered log-volatilities
#' (i.e., the log-vols minus the unconditional means \code{dhs_mean})
#' @param h_phi the \code{p x 1} vector of previous AR(1) coefficient(s)
#' @param h_sigma_eta_t the \code{T x p} matrix of log-vol innovation standard deviations
#' @param prior_dhs_phi the parameters of the prior for the log-volatilty AR(1) coefficient \code{dhs_phi};
#' either \code{NULL} for uniform on [-1,1] or a 2-dimensional vector of (shape1, shape2) for a Beta prior
#' on \code{[(dhs_phi + 1)/2]}
#'
#' @return \code{p x 1} vector of sampled AR(1) coefficient(s)
#'
#' @note For the standard AR(1) case, \code{p = 1}. However, the function applies more
#' generally for sampling \code{p > 1} independent AR(1) processes (jointly).
#'
#' @import truncdist
#' @export
#' 
#----------------------------------------------------------------------------
sampleAR1 = function(h_yc, h_phi, h_sigma_eta_t, prior_dhs_phi = NULL){
  
  # Compute dimensions:
  n = nrow(h_yc); p = ncol(h_yc)
  
  # Loop over the j=1:p
  for(j in 1:p){
    
    # Compute "regression" terms for dhs_phi_j:
    y_ar = h_yc[-1,j]/h_sigma_eta_t[,j] # Standardized "response"
    x_ar = h_yc[-n,j]/h_sigma_eta_t[,j] # Standardized "predictor"
    
    # Using Beta distribution:
    if(!is.null(prior_dhs_phi)){
      
      # Check to make sure the prior params make sense
      if(length(prior_dhs_phi) != 2) stop('prior_dhs_phi must be a numeric vector of length 2')
      
      dhs_phi01 = (h_phi[j] + 1)/2 # ~ Beta(prior_dhs_phi[1], prior_dhs_phi[2])
      
      # Slice sampler when using Beta prior:
      dhs_phi01 = uni.slice(dhs_phi01, g = function(x) {
        -0.5*sum((y_ar - (2*x - 1)*x_ar)^2) +
          dbeta(x, shape1 = prior_dhs_phi[1], shape2 = prior_dhs_phi[2], log = TRUE)
      }, lower = 0, upper = 1)[1]#}, lower = 0.005, upper = 0.995)[1] #
      
      h_phi[j] = 2*dhs_phi01 - 1
      
    } else {
      # For h_phi ~ Unif(-1, 1), the posterior is truncated normal
      h_phi[j] = rtrunc(n = 1, spec = 'norm',
                        a = -1, b = 1,
                        mean = sum(y_ar*x_ar)/sum(x_ar^2),
                        sd = 1/sqrt(sum(x_ar^2)))
    }
  }
  h_phi
}
#----------------------------------------------------------------------------
#' Sample the dynamic shrinkage process parameters
#'
#' Compute one draw for each of the parameters in the dynamic shrinkage process
#' for the special case in which the shrinkage parameter \code{kappa ~ Beta(alpha, beta)}
#' with \code{alpha = beta}. The primary example is the dynamic horseshoe process with
#' \code{alpha = beta = 1/2}.
#'
#' @param omega \code{T x p} matrix of evolution errors
#' @param evolParams list of parameters to be updated (see Value below)
#' @param sigma_e the observation error standard deviation; for (optional) scaling purposes
#' @param prior_dhs_phi the parameters of the prior for the log-volatilty AR(1) coefficient \code{dhs_phi};
#' either \code{NULL} for uniform on [-1,1] or a 2-dimensional vector of (shape1, shape2) for a Beta prior
#' on \code{[(dhs_phi + 1)/2]}
#' @param alphaPlusBeta For the symmetric prior kappa ~ Beta(alpha, beta) with alpha=beta,
#' specify the sum [alpha + beta]
#' @return List of relevant components:
#' \itemize{
#' \item the \code{T x p} evolution error standard deviations \code{sigma_wt},
#' \item the \code{T x p} log-volatility \code{ht}, the \code{p x 1} log-vol unconditional mean(s) \code{dhs_mean},
#' \item the \code{p x 1} log-vol AR(1) coefficient(s) \code{dhs_phi},
#' \item the \code{T x p} log-vol innovation standard deviations \code{sigma_eta_t} from the Polya-Gamma priors,
#' \item the \code{p x 1} initial log-vol SD \code{sigma_eta_0},
#' \item and the mean of log-vol means \code{dhs_mean0} (relevant when \code{p > 1})
#' }
#'
#' @note The priors induced by \code{prior_dhs_phi} all imply a stationary (log-) volatility process.
#'
#' @import BayesLogit
#' @export
sampleDSP = function(omega, evolParams, sigma_e = 1, prior_dhs_phi = c(10,2),   alphaPlusBeta = 1){
  
  # Store the DSP parameters locally:
  ht = evolParams$ht; dhs_mean = evolParams$dhs_mean; 
  dhs_phi = evolParams$dhs_phi; sigma_eta_t = evolParams$sigma_eta_t; 
  sigma_eta_0 = evolParams$sigma_eta_0; dhs_mean0 = evolParams$dhs_mean0
  
  # "Local" number of time points
  ht = as.matrix(ht)
  n = nrow(ht); p = ncol(ht)
  
  # Sample the log-volatilities using AWOL sampler
  ht = sampleLogVols(h_y = omega, h_prev = ht, h_mu = dhs_mean, h_phi=dhs_phi, h_sigma_eta_t = sigma_eta_t, h_sigma_eta_0 = sigma_eta_0)
  
  # Compute centered log-vols for the samplers below:
  ht_tilde = ht - tcrossprod(rep(1,n), dhs_mean)
  
  # Sample AR(1) parameters
  # Note: dhs_phi = 0 means non-dynamic HS, while dhs_phi = 1 means RW, in which case we don't sample either
  if(!all(dhs_phi == 0) && !all(dhs_phi == 1)) dhs_phi = sampleAR1(h_yc = ht_tilde, h_phi = dhs_phi, h_sigma_eta_t = sigma_eta_t, prior_dhs_phi = prior_dhs_phi)
  
  # Sample the evolution error SD of log-vol (i.e., Polya-Gamma mixing weights)
  eta_t = ht_tilde[-1,] - tcrossprod(rep(1,n-1), dhs_phi)*ht_tilde[-n, ]       # Residuals
  sigma_eta_t = matrix(1/sqrt(rpg(num = (n-1)*p, h = alphaPlusBeta, z = eta_t)), nc = p) # Sample
  sigma_eta_0 = 1/sqrt(rpg(num = p, h = 1, z = ht_tilde[1,]))                # Sample the inital
  
  # Sample the unconditional mean(s), unless dhs_phi = 1 (not defined)
  if(!all(dhs_phi == 1)){
    if(p > 1){
      # Assume a hierarchy of the global shrinkage params across j=1,...,p
      muSample = sampleLogVolMu(h = ht, h_mu = dhs_mean, h_phi = dhs_phi, h_sigma_eta_t = sigma_eta_t, h_sigma_eta_0 = sigma_eta_0, h_log_scale = dhs_mean0);
      dhs_mean = mlLuSample$dhs_mean
      dhs_mean0 = sampleLogVolMu0(h_mu = dhs_mean, h_mu0 = dhs_mean0, dhs_mean_prec_j = muSample$dhs_mean_prec_j, h_log_scale = log(sigma_e^2))
    } else {
      # p = 1
      muSample = sampleLogVolMu(h = ht, h_mu = dhs_mean, h_phi = dhs_phi, h_sigma_eta_t = sigma_eta_t, h_sigma_eta_0 = sigma_eta_0, h_log_scale = log(sigma_e^2));
      dhs_mean = dhs_mean0 = muSample$dhs_mean # save dhs_mean0 = dhs_mean for coding convenience later
    }
  } else {dhs_mean = rep(0, p); dhs_mean0 = 0} # When RW for log-vols, fix unconditional mean for identifiability
  
  # Evolution error SD:
  sigma_wt = exp(ht/2)
  
  # Return the same list, but with the new values
  list(sigma_wt = sigma_wt, ht = ht, dhs_mean = dhs_mean, dhs_phi = dhs_phi, sigma_eta_t = sigma_eta_t, sigma_eta_0 = sigma_eta_0, dhs_mean0 = dhs_mean0)
}
#----------------------------------------------------------------------------
#' Sample the AR(1) unconditional means
#'
#' Compute one draw of the unconditional means in an AR(1) model with Gaussian innovations
#' and time-dependent innovation variances. In particular, we use the sampler for the
#' log-volatility AR(1) process with the parameter-expanded Polya-Gamma sampler. The sampler also applies
#' to a multivariate case with independent components.
#'
#' @param h the \code{T x p} matrix of log-volatilities
#' @param h_mu the \code{p x 1} vector of previous means
#' @param h_phi the \code{p x 1} vector of AR(1) coefficient(s)
#' @param h_sigma_eta_t the \code{T x p} matrix of log-vol innovation standard deviations
#' @param h_sigma_eta_0 the standard deviations of initial log-vols
#' @param h_log_scale prior mean from scale mixture of Gaussian (Polya-Gamma) prior, e.g. log(sigma_e^2) or dhs_mean0
#'
#' @return a list containing
#' \itemize{
#' \item the sampled mean(s) \code{dhs_mean} and
#' \item the sampled precision(s) \code{dhs_mean_prec_j} from the Polya-Gamma parameter expansion
#'}
#'
#' @import BayesLogit
#' @export
sampleLogVolMu = function(h, h_mu, h_phi, h_sigma_eta_t, h_sigma_eta_0, h_log_scale = 0){
  
  # Compute "local" dimensions:
  n = nrow(h); p = ncol(h)
  
  # Sample the precision term(s)
  dhs_mean_prec_j = rpg(num = p, h = 1, z = h_mu - h_log_scale)
  
  # Now, form the "y" and "x" terms in the (auto)regression
  y_mu = (h[-1,] - tcrossprod(rep(1,n-1), h_phi)*h[-n,])/h_sigma_eta_t;
  x_mu = tcrossprod(rep(1,n-1), 1 - h_phi)/h_sigma_eta_t
  
  # Include the initial sd?
  #if(!is.null(h_sigma_eta_0)){y_mu = rbind(h[1,]/h_sigma_eta_0, y_mu); x_mu = rbind(1/h_sigma_eta_0, x_mu)}
  y_mu = rbind(h[1,]/h_sigma_eta_0, y_mu);
  x_mu = rbind(1/h_sigma_eta_0, x_mu)
  
  # Posterior SD and mean:
  postSD = 1/sqrt(colSums(x_mu^2) + dhs_mean_prec_j)
  postMean = (colSums(x_mu*y_mu) + h_log_scale*dhs_mean_prec_j)*postSD^2
  dhs_mean = rnorm(n = p, mean = postMean, sd = postSD)
  
  list(dhs_mean = dhs_mean, dhs_mean_prec_j = dhs_mean_prec_j)
}
#----------------------------------------------------------------------------
#' Sample the latent log-volatilities
#'
#' Compute one draw of the log-volatilities using a discrete mixture of Gaussians
#' approximation to the likelihood (see Omori, Chib, Shephard, and Nakajima, 2007)
#' where the log-vols are assumed to follow an AR(1) model with time-dependent
#' innovation variances. More generally, the code operates for \code{p} independent
#' AR(1) log-vol processes to produce an efficient joint sampler in \code{O(Tp)} time.
#'
#' @param h_y the \code{T x p} matrix of data, which follow independent SV models
#' @param h_prev the \code{T x p} matrix of the previous log-vols
#' @param h_mu the \code{p x 1} vector of log-vol unconditional means
#' @param h_phi the \code{p x 1} vector of log-vol AR(1) coefficients
#' @param h_sigma_eta_t the \code{T x p} matrix of log-vol innovation standard deviations
#' @param h_sigma_eta_0 the \code{p x 1} vector of initial log-vol innovation standard deviations
#'
#' @return \code{T x p} matrix of simulated log-vols
#'
#' @note For Bayesian trend filtering, \code{p = 1}. More generally, the sampler allows for
#' \code{p > 1} but assumes (contemporaneous) independence across the log-vols for \code{j = 1,...,p}.
#'
#' @import Matrix
#' @import BayesLogit
sampleLogVols = function(h_y, h_prev, h_mu, h_phi, h_sigma_eta_t, h_sigma_eta_0){
  
  # Compute dimensions:
  h_prev = as.matrix(h_prev) # Just to be sure (T x p)
  n = nrow(h_prev); p = ncol(h_prev)
  
  # Mixture params: mean, variance, and weights
  # Kim, Shephard, Chib (1998) 7-component mixture:
  #m_st  = c(-11.40039, -5.24321, -9.83726, 1.50746,  -0.65098, 0.52478,  -2.35859)
  #v_st2 = c(5.795960,  2.613690, 5.179500, 0.167350, 0.640090, 0.340230, 1.262610)
  #q     = c(0.007300,  0.105560, 0.000020, 0.043950, 0.340010, 0.245660, 0.257500)
  
  # Omori, Chib, Shephard, Nakajima (2007) 10-component mixture:
  m_st  = c(1.92677, 1.34744, 0.73504, 0.02266, -0.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000)
  v_st2 = c(0.11265, 0.17788, 0.26768, 0.40611,  0.62699,  0.98583,  1.57469,  2.54498,  4.16591,   7.33342)
  q     = c(0.00609, 0.04775, 0.13057, 0.20674,  0.22715,  0.18842,  0.12047,  0.05591,  0.01575,   0.00115)
  
  # Add an offset: common for all times, but distict for each j=1,...,p
  yoffset = tcrossprod(rep(1,n),
                       apply(as.matrix(h_y), 2,
                             function(x) any(x^2 < 10^-16)*max(10^-8, mad(x)/10^6)))
  
  # This is the response in our DLM, log(y^2)
  ystar = log(h_y^2 + yoffset)
  
  # Sample the mixture components
  #z = draw.indicators(res = ystar-h_prev, nmix = list(m = m_st, v = v_st2, p = q))
  z = sapply(ystar-h_prev, ncind, m_st, sqrt(v_st2), q)
  
  # Subset mean and variances to the sampled mixture components; (n x p) matrices
  m_st_all = matrix(m_st[z], nr=n); v_st2_all = matrix(v_st2[z], nr=n)
  
  # Joint AWOL sampler for j=1,...,p:
  
  # Constant (but j-specific) mean
  h_mu_all = tcrossprod(rep(1,n), h_mu)
  
  # Constant (but j-specific) AR(1) coef
  h_phi_all = tcrossprod(rep(1,n), h_phi)
  
  # Linear term:
  linht = matrix((ystar - m_st_all - h_mu_all)/v_st2_all)
  
  # Evolution precision matrix (n x p)
  evol_prec_mat = matrix(0, nr = n, nc = p);
  evol_prec_mat[1,] = 1/h_sigma_eta_0^2;
  evol_prec_mat[-1,] = 1/h_sigma_eta_t^2;
  
  # Lagged version, with zeros as appropriate (needed below)
  evol_prec_lag_mat = matrix(0, nr = n, nc = p);
  evol_prec_lag_mat[1:(n-1),] = evol_prec_mat[-1,]
  
  # Diagonal of quadratic term:
  Q_diag = matrix(1/v_st2_all +  evol_prec_mat + h_phi_all^2*evol_prec_lag_mat)
  
  # Off-diagonal of quadratic term:
  Q_off = matrix(-h_phi_all*evol_prec_lag_mat)[-(n*p)]
  
  # Quadratic term:
  QHt_Matrix = bandSparse(n*p, k = c(0,1), diag = list(Q_diag, Q_off), symm = TRUE)
  #QHt_Matrix = as.spam.dgCMatrix(as(bandSparse(n*p, k = c(0,1), diag = list(Q_diag, Q_off), symm = TRUE),"dgCMatrix"))
  
  # Cholesky:
  chQht_Matrix = Matrix::chol(QHt_Matrix)
  
  # Sample the log-vols:
  hsamp = h_mu_all + matrix(Matrix::solve(chQht_Matrix,Matrix::solve(Matrix::t(chQht_Matrix), linht) + rnorm(length(linht))), nr = n)
  #hsamp = h_mu_all +matrix(rmvnorm.canonical(n = 1, b = linht, Q = QHt_Matrix, Rstruct = cholDSP0))
  
  
  # Return the (uncentered) log-vols
  hsamp
}
#----------------------------------------------------------------------------
#' Univariate Slice Sampler from Neal (2008)
#'
#' Compute a draw from a univariate distribution using the code provided by
#' Radford M. Neal. The documentation below is also reproduced from Neal (2008).
#'
#' @param x0    Initial point
#' @param g     Function returning the log of the probability density (plus constant)
#' @param w     Size of the steps for creating interval (default 1)
#' @param m     Limit on steps (default infinite)
#' @param lower Lower bound on support of the distribution (default -Inf)
#' @param upper Upper bound on support of the distribution (default +Inf)
#' @param gx0   Value of g(x0), if known (default is not known)
#'
#' @return  The point sampled, with its log density attached as an attribute.
#'
#' @note The log density function may return -Inf for points outside the support
#' of the distribution.  If a lower and/or upper bound is specified for the
#' support, the log density function will not be called outside such limits.
#' 
uni.slice = function (x0, g, w=1, m=Inf, lower=-Inf, upper=+Inf, gx0=NULL)
{
  # Check the validity of the arguments.
  
  if (!is.numeric(x0) || length(x0)!=1
      || !is.function(g)
      || !is.numeric(w) || length(w)!=1 || w<=0
      || !is.numeric(m) || !is.infinite(m) && (m<=0 || m>1e9 || floor(m)!=m)
      || !is.numeric(lower) || length(lower)!=1 || x0<lower
      || !is.numeric(upper) || length(upper)!=1 || x0>upper
      || upper<=lower
      || !is.null(gx0) && (!is.numeric(gx0) || length(gx0)!=1))
  {
    stop ("Invalid slice sampling argument")
  }
  
  # Keep track of the number of calls made to this function.
  #uni.slice.calls <<- uni.slice.calls + 1
  
  # Find the log density at the initial point, if not already known.
  
  if (is.null(gx0))
  { #uni.slice.evals <<- uni.slice.evals + 1
    gx0 <- g(x0)
  }
  
  # Determine the slice level, in log terms.
  
  logy <- gx0 - rexp(1)
  
  # Find the initial interval to sample from.
  
  u <- runif(1,0,w)
  L <- x0 - u
  R <- x0 + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff
  
  # Expand the interval until its ends are outside the slice, or until
  # the limit on steps is reached.
  
  if (is.infinite(m))  # no limit on number of steps
  {
    repeat
    { if (L<=lower) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
    }
    
    repeat
    { if (R>=upper) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
    }
  }
  
  else if (m>1)  # limit on steps, bigger than one
  {
    J <- floor(runif(1,0,m))
    K <- (m-1) - J
    
    while (J>0)
    { if (L<=lower) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
      J <- J - 1
    }
    
    while (K>0)
    { if (R>=upper) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
      K <- K - 1
    }
  }
  
  # Shrink interval to lower and upper bounds.
  
  if (L<lower)
  { L <- lower
  }
  if (R>upper)
  { R <- upper
  }
  
  # Sample from the interval, shrinking it on each rejection.
  
  repeat
  {
    x1 <- runif(1,L,R)
    
    #uni.slice.evals <<- uni.slice.evals + 1
    gx1 <- g(x1)
    
    if (gx1>=logy) break
    
    if (x1>x0)
    { R <- x1
    }
    else
    { L <- x1
    }
  }
  
  # Return the point sampled, with its log density attached as an attribute.
  
  attr(x1,"log.density") <- gx1
  return (x1)
  
}


