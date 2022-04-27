# wrappers to unify probability estimation interface and output ----------------

sov = function(mu, Sigma, lb, ub, n_batch_mc) {
  gb <- tlrmvnmvt::GenzBretz(n_batch_mc)
  
  prob <- tlrmvnmvt::pmvn(lb, ub, mu, Sigma, algorithm = gb)
  result <- data.frame(variable = c("estimate", "error"), 
                       value = c(prob, attr(prob, "error")))
  return(result)
}

# block_size set to sqrt(d), as in Cao et al 2021
tlr <- function(mu, Sigma, lb, ub, n_batch_mc, block_size = NULL,
                epsl = 1e-4) {
  d <- length(mu)
  if (is.null(block_size))
    m = max(floor(sqrt(d)), 4)
  else 
    m = block_size
  tlr <- tlrmvnmvt::TLRQMC(N = n_batch_mc, m = m, epsl = epsl)
  prob <- tlrmvnmvt::pmvn(lb, ub, mu, Sigma, algorithm = tlr)
  result <- 
    data.frame(variable = c("estimate", "error"),
               value = c(prob, attr(prob, "error")))
  return(result)
}

epmgp <- function(mu, Sigma, lb, ub) {
  prob <- epmgpr::pmvn(lb, ub, mu, Sigma)
  return(data.frame(variable = "estimate", value = prob))
}

ghk <- function(mu, Sigma, lb, ub, n_batch_mc, n_est = 10) {
  
  ests = rep(NA, n_est)
  for (i in 1:n_est) {
    result <- gcKrig::mvnintGHK(mu, Sigma, lb, ub, nrep = n_batch_mc, 
                                reorder = TRUE, log = FALSE)  
    ests[i] = result$value
  }
  return(data.frame(
    variable = c("estimate", "error"),
    value = c(mean(ests), sd(ests) / sqrt(n_est))
  ))
}

# from the paper: 
# n_sub_samples = 16
# n_sub_skip = 10
lcg <- function(mu, Sigma, lb, ub, 
                n_sub_samples, n_hdr_samples, 
                n_sub_skip = 1, n_hdr_skip = 1,
                domain_fraction = .5, 
                n_est = 10) {
  prob <- lincongauss::pmvn(mu, Sigma, lb, ub, 
                            n_sub_samples = n_sub_samples,
                            n_sub_skip = n_sub_skip,
                            n_hdr_samples = n_hdr_samples,
                            n_hdr_skip = n_hdr_skip,
                            domain_fraction = domain_fraction, 
                            n_est = n_est)
  return(data.frame(
    variable = c("estimate", "error"), 
    value = c(prob, attr(prob, "error"))
  ))
}

met = function(mu, Sigma, lb, ub, n_batch_mc, n_est = 10) {
  prob = met::pmvn(mu, Sigma, lb, ub, n = n_batch_mc, n_est = n_est)
  return(data.frame(
    variable = c("estimate", "relerror", "upvnd"),
    value = c(prob, attr(prob, "relErr"), attr(prob, "upbnd"))
  ))
}

bvcdn = function(mu, Sigma, lb, ub, tol = 1e-4) {
  problem_d = length(mu)
  vorder = hccmvn::recurUniandBlkcmb(
    covM = Sigma, a = lb - mu, b = ub - mu, bsz_ = 1)
  Sigma = Sigma[vorder, vorder]
  lb = lb[vorder]
  ub = ub[vorder]
  mu = mu[vorder]
  prob = hccmvn::hccmvn(Sigma, lb - mu, ub - mu, m=problem_d, d=2, tol=tol)
  return(data.frame(
    variable = "estimate",
    value = prob
  ))
}

hblkcdn = function(mu, Sigma, lb, ub, cond_size, block_size, tol = 1e-4) {
  problem_d = length(mu)
  vorder = hccmvn::recurUniandBlkcmb(
    covM = Sigma, a = lb - mu, b = ub - mu, bsz_ = block_size)
  Sigma = Sigma[vorder, vorder]
  lb = lb[vorder]
  ub = ub[vorder]
  mu = mu[vorder]
  prob = hccmvn::hccmvn(Sigma, lb - mu, ub - mu, m=block_size, d=cond_size,
                        tol=tol)
  return(data.frame(
    variable = "estimate",
    value = prob
  ))
}


