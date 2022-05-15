# wrappers to unify sampling function interfaces

gibbs_cov <- function(n, mu, Sigma, lb, ub, initial, burnin = 0) {
  tmvtnorm::rtmvnorm(n, mu, Sigma, lb, ub, algorithm = "gibbs",
                     burn.in.samples = burnin, start.value = initial)
}

gibbs_prec = function(n, mu, Sigma, lb, ub, initial, burnin = 0) {
  Prec = chol2inv(chol(Sigma))
  tmvtnorm::rtmvnorm(n,
                     mu, H = Prec, 
                     lower = lb, upper = ub, 
                     algorithm = "gibbs", 
                     burn.in.samples = burnin, 
                     start.value = initial)
}

gibbs_ry2004 = function(n, mu, Sigma, lb, ub, A = diag(length(mu)), 
                        initial, burnin = 0, tuvn_sampler) {
  sampletmvn::sample_gibbs_ry2004(n, mu, Sigma, lb, ub, A,
                                  init = initial, 
                                  burn =  burnin, 
                                  thin = 0, 
                                  tuvn_sampler = tuvn_sampler)
}

gibbs_lg2015 = function(n, mu, Sigma, lb, ub, A = diag(length(mu)),
                        initial, burnin = 0, tuvn_sampler) {
  sampletmvn::sample_gibbs_lg2015(n, mu, Sigma, lb, ub, A, 
                                  init = initial,
                                  burn = burnin,
                                  thin = 0,
                                  tuvn_sampler = tuvn_sampler)  
}

ehmc <- function(n, mu, Sigma, lb, ub, A = diag(length(mu)), initial, burnin = 0) {
  d <- length(mu)
  Prec <- chol2inv(chol(Sigma))
  r <- as.vector(Prec %*% mu)
  inf_idx <- c(is.infinite(lb), is.infinite(ub))
  f <- rbind(A, -A)[!inf_idx, , drop = FALSE]
  g <- c(lb, ub)[!inf_idx]
  tmg::rtmg(n, Prec, r, initial, f, g, burn.in = burnin)
}

rsm = function(n, mu, Sigma, lb, ub, A = NULL, initial = NULL) {
  sampletmvn::sample_rsm(n, mu, Sigma, lb, ub, A)
}

met = function(n, mu, Sigma, lb, ub, A = NULL, initial = NULL) {
  met::rtmvn(n, mu, Sigma, lb, ub, A)
}

met_mp = function(n, mu, Sigma, lb, ub, A, initial = NULL) {
  
  G = MASS::ginv(A)
  m = nrow(A)
  Amu = A %*% mu
  lb_t = lb - Amu
  ub_t = ub - Amu
  Sigma_t = A %*% Sigma %*% t(A)
  mu_t = rep(0, m)
  
  z = met::rtmvn(n, mu_t, Sigma_t, lb_t, ub_t)
  
  x = G %*% t(z + matrix(rep(Amu, n), nrow = n, byrow = TRUE))
  
  return(t(x))
}

epess = function(n, mu, Sigma, lb, ub, A = NULL, initial, 
                 n_perthresh, n_slices, burnin = 0) {
  epmgpr::rtmvn(n, mu, Sigma, lb, ub, A, 
                initial = initial, 
                n_perthresh = n_perthresh, n_slices = n_slices, burnin = burnin)
}

liness = function(n, mu, Sigma, lb, ub, A = NULL, initial) {
  lincongauss::rtmvn(n, mu, Sigma, lb, ub, A, 
                     x_init = initial)
}

rhmc = function(n, mu, Sigma, lb, ub, initial, traj_length, burnin = 0) {
  sampletmvn::sample_rhmc(n, mu, Sigma, lb, ub, initial = initial, burnin = burnin,
                          traj_length = traj_length)
}

ghk = function(n, mu, Sigma, lb, ub, initial = NULL) {
  d = length(mu)
  result = tmvnsim::tmvnsim(n, d, lower = lb, upper = ub, 
                            means = mu, sigma = Sigma)
  samples <- result$samp
  weights <- result$wts
  idx <- sample(1:nrow(samples), size = n, replace = TRUE, prob = weights)
  return(samples[idx, ])
}







