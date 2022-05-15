library(glasso)

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

probit_mcmc_glasso = function(X, Y, sampler, sampler_params, penalty = .01,
                              Sigma_init = NULL, beta_init = NULL, 
                              tol = 1e-5, max_iter = 200, n_mc_samples = 100,
                              corr = FALSE, pen_diag = FALSE,
                              verbose = FALSE) {
  
  n_choices = ncol(Y)
  n_obs = nrow(Y)
  p = ncol(X[[1]])
  
  lb = rep(0, n_choices-1)
  ub = rep(Inf, n_choices-1)
  
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
  dPrec = Inf
  glasso_llik = rep(NA, max_iter)
  
  constraints = apply(Y, 1, max_polytope, simplify = FALSE)
  initial_points = apply(Y, 1, function(y) 2*y - 1, simplify = FALSE)
  
  Sigma_history = vector(mode = "list", length = max_iter)
  Prec_history = vector(mode = "list", length = max_iter)
  
  # while (dSigma > tol & iter <= max_iter) {
  while (dPrec > tol & iter <= max_iter) {
    
    Prec = chol2inv(chol(Sigma))
    E_sample_cov = matrix(0, nrow = n_choices, ncol = n_choices)
    
    # E step
    for (i in 1:n_obs) {
      
      mu = X[[i]] %*% beta
      A = constraints[[i]]
      initial = initial_points[[i]]
      
      problem_params = list(n = n_mc_samples, 
                            mu = mu, Sigma = Sigma, 
                            lb = lb, ub = ub, A = A,
                            initial = initial)
      params = c(problem_params, sampler_params)
      
      samples = do.call(sampler, params)
      
      mc_est = matrix(0, nrow = n_choices, ncol = n_choices)
      for (j in 1:n_mc_samples) {
        mc_est = mc_est + tcrossprod(samples[j, ] - mu)
      }
      mc_est = mc_est / n_mc_samples
      
      E_sample_cov = E_sample_cov + mc_est
    }
    
    E_sample_cov = E_sample_cov / n_obs
    E_sample_cor = cov2cor(E_sample_cov)
    
    # M step
    glasso_result = glasso(E_sample_cor, penalty, nobs = n_obs, 
                           penalize.diagonal = pen_diag,
                           start = "warm", w.init = Sigma, wi.init = Prec)
    glasso_llik[iter] = glasso_result$loglik
    
    if (verbose)
      print(paste0("Iter ", iter, ", llik: ", glasso_llik[iter]))
    
    Sigma_new = glasso_result$w
    Prec_new = glasso_result$wi
    
    if (corr) {
      Sigma_new = cov2cor(Sigma_new)
      Prec_new = cov2cor(Prec_new)
    }
    
    Sigma_history[[iter]] = Sigma_new
    Prec_history[[iter]] = Prec_new
    
    dSigma = max(abs(Sigma_new  - Sigma))
    dPrec = max(abs(Prec_new - Prec))
    iter = iter + 1
    Sigma = Sigma_new
    Prec = Prec_new
  }
  
  # keep only elements before convergence
  glasso_llik = glasso_llik[1:iter]
  Sigma_history = Sigma_history[1:iter]
  Prec_history = Prec_history[1:iter]
  
  return(list(
    Sigma = Sigma,
    Prec = Prec,
    beta = beta,
    iters = iter,
    gllaso_llik = glasso_llik,
    algorithm = algorithm,
    n_mc_samples = n_mc_samples,
    Sigma_history = Sigma_history,
    Prec_history = Prec_history
  ))
}