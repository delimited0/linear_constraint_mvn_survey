library(CovTools)

max_polytope = function(y) {
  n = length(y)
  idx = which(y == 1)
  
  D = -1*diag(n-1)
  if (idx == 1)
    A = cbind(1, D)
  else if (idx == n)
    A = cbind(D, 1)
  else
    A = cbind(D[, 1:(idx-1)], rep(1, n-1), D[, idx:(n-1)])
  
  return(A)
}


probit_sparse_cov = function(X, Y, sampler, initial,
                             penalty = .01,
                             Sigma_init = NULL, beta_init = NULL, 
                             tol = 1e-5, max_iter = 200, n_mc_samples = 100,
                             corr = FALSE, pen_diag = FALSE) {
  
  n_choices = ncol(Y)
  n_obs = nrow(Y)
  p = ncol(X[[1]])
  
  lb = rep(0, n_choices)
  ub = rep(Inf, n_choices)
  
  if (is.null(Sigma_init))
    Sigma = diag(n_choices)
  else 
    Sigma = Sigma_init
  
  if (is.null(beta_init))
    beta = rep(0, p)
  else
    beta = beta_init
  
  iter = 1
  dSigma = Inf
  glasso_llik = rep(NA, max_iter)
  
  chain_init = initial
  constraints = apply(Y, 1, max_polytope, simplify = FALSE)
  
  E_cov = 
  
  while (dSigma > tol & iter <= max_iter) {
    
    for (i in 1:n_obs) {
      
      mu = X[[i]] %*% beta
      A = constraints[[i]]
      
      samples = sampler(n_mc_samples, mu, Sigma, lb, ub, A, 
                        initial = chain_init)
      
      mc_est = matrix(0, nrow = n_choices, ncol = n_choices)
      for (j in 1:n_mc_samples) {
        mc_est = mc_est + tcrossprod(samples[j, ] - mu)
      }
      mc_est = mc_est / mc_samples
      
      E_sample_cov = E_sample_cov + mc_est
    }
    
    E_sample_cov = E_sample_cov / n_obs
    
    
    
  }
}