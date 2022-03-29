# wrappers to unify sampling function interfaces

gibbs_cov <- function(n, mu, Sigma, lb, ub, initial, burnin) {
  tmvtnorm::rtmvnorm(n, mu, Sigma, lb, ub, algorithm = "gibbs",
                     burn.in.samples = burnin, start.value = initial)
}

gibbs_prec = function(n, mu, Sigma, lb, ub, initial, burnin) {
  Prec = chol2inv(chol(Sigma))
  tmvtnorm::rtmvnorm(n,
                     mu, H = Prec, 
                     lower = lb, upper = ub, 
                     algorithm = "gibbs", 
                     burn.in.samples = burnin, 
                     start.value = initial)
}

gibbs_ry2004 = function(n, mu, Sigma, lb, ub, A, initial, burnin, tuvn_sampler) {
  sampletmvn::sample_gibbs_ry2004(n, mu, Sigma, A, lb, ub, 
                                  init = initial, 
                                  burn =  burnin, 
                                  thin = 0, 
                                  tuvn_sampler = tuvn_sampler)
}

ehmc <- function(n, mu, Sigma, lb, ub, A, initial, burnin) {
  d <- length(mu)
  Prec <- chol2inv(chol(Sigma))
  r <- as.vector(Prec %*% mu)
  inf_idx <- c(is.infinite(lb), is.infinite(ub))
  f <- rbind(diag(d), -diag(d))[!inf_idx, , drop = FALSE]
  g <- c(-lb, ub)[!inf_idx]
  tmg::rtmg(n, Prec, r, initial, f, g, burn.in = burnin)
}



