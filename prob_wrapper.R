# wrappers to unify probability estimation interface and output ----------------


# Genz and Bretz QRSVN ----------------------------------------------------
sov = function(mu, Sigma, lb, ub, n_batch_mc) {
  gb <- tlrmvnmvt::GenzBretz(n_batch_mc)
  
  prob <- tlrmvnmvt::pmvn(lb, ub, mu, Sigma, algorithm = gb)
  return(prob)
}


# Tile low rank approximation ---------------------------------------------
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
  return(prob)
}


# Expectation propagation -------------------------------------------------
epmgp <- function(mu, Sigma, lb, ub) {
  prob <- epmgpr::pmvn(lb, ub, mu, Sigma)
  return(prob)
}


# GHK ---------------------------------------------------------------------
ghk <- function(mu, Sigma, lb, ub, n_batch_mc, n_est = 10) {
  
  ests = rep(NA, n_est)
  for (i in 1:n_est) {
    result <- gcKrig::mvnintGHK(mu, Sigma, lb, ub, nrep = n_batch_mc, 
                                reorder = TRUE, log = FALSE)  
    ests[i] = result$value
  }
  
  prob = mean(ests)
  attr(prob, "error") = sd(ests) / sqrt(n_est)
  
  return(prob)
}



# LCG ---------------------------------------------------------------------
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
  return(prob)
}


# Minimax exponential tilting ---------------------------------------------
met = function(mu, Sigma, lb, ub, n_batch_mc, n_est = 10) {
  prob = met::pmvn(mu, Sigma, lb, ub, n = n_batch_mc, n_est = n_est)
  return(prob)
}



# univariate conditioning -------------------------------------------------
uvcdn = function(mu, Sigma, lb, ub, tol = 1e-4) {
  
  problem_d = length(mu)
  odr = hccmvn::recurUniandBlkcmb(Sigma, lb - mu, ub - mu, problem_d)+1
  
  # replace infinity with large numbers
  sds = sqrt(diag(Sigma))
  lb_inf_idx = is.infinite(lb)
  ub_inf_idx = is.infinite(ub)
  
  lb[lb_inf_idx] = -10 * sds[lb_inf_idx]
  ub[ub_inf_idx] = 10 * sds[ub_inf_idx]
  
  prob = hccmvn::hccmvn(
    covM = Sigma[odr, odr],
    a = (lb - mu)[odr],
    b = (ub - mu)[odr],
    m = problem_d, 
    d = 1, 
    tol = 1e-4)
  
  return(prob)
}

# bivariate conditioning -------------------------------------------------
bvcdn = function(mu, Sigma, lb, ub, tol = 1e-4) {
  
  problem_d = length(mu)
  odr = hccmvn::recurUniandBlkcmb(Sigma, lb - mu, ub - mu, problem_d)+1
  
  # replace infinity with large numbers
  sds = sqrt(diag(Sigma))
  lb_inf_idx = is.infinite(lb)
  ub_inf_idx = is.infinite(ub)
  
  lb[lb_inf_idx] = -10 * sds[lb_inf_idx]
  ub[ub_inf_idx] = 10 * sds[ub_inf_idx]
  
  prob = hccmvn::hccmvn(
    covM = Sigma[odr, odr],
    a = (lb - mu)[odr],
    b = (ub - mu)[odr],
    m = problem_d, 
    d = 2, 
    tol = 1e-4)
  
  return(prob)
}

# d-variate conditioning --------------------------------------------------
dvcdn = function(mu, Sigma, lb, ub, cond_size, tol = 1e-4) {
  
  problem_d = length(mu)
  odr = hccmvn::recurUniandBlkcmb(Sigma, lb - mu, ub - mu, problem_d)+1
  
  # replace infinity with large numbers
  sds = sqrt(diag(Sigma))
  lb_inf_idx = is.infinite(lb)
  ub_inf_idx = is.infinite(ub)
  
  lb[lb_inf_idx] = -10 * sds[lb_inf_idx]
  ub[ub_inf_idx] = 10 * sds[ub_inf_idx]
  
  prob = hccmvn::hccmvn(
    covM = Sigma[odr, odr],
    a = (lb - mu)[odr],
    b = (ub - mu)[odr],
    m = problem_d, 
    d = cond_size, 
    tol = 1e-4)
  
  return(prob)
}

# Hierarchical block conditioning -----------------------------------------

# problem dimension should be equal to block size * 2^{int}
hblkcdn = function(mu, Sigma, lb, ub, cond_size, block_size, tol = 1e-4) {
  problem_d = length(mu)
  odr = hccmvn::recurUniandBlkcmb(Sigma, lb - mu, ub - mu, block_size)+1
  
  # replace infinity with large numbers
  sds = sqrt(diag(Sigma))
  lb_inf_idx = is.infinite(lb)
  lb_inf_sign = sign(lb[lb_inf_idx])
  ub_inf_idx = is.infinite(ub)
  ub_inf_sign = sign(ub[ub_inf_idx])
  
  lb[lb_inf_idx] = lb_inf_sign * 10 * sds[lb_inf_idx]
  ub[ub_inf_idx] = ub_inf_sign * 10 * sds[ub_inf_idx]
  
  prob = hccmvn::hccmvn(
    covM = Sigma[odr, odr],
    a = (lb - mu)[odr],
    b = (ub - mu)[odr], 
    m = block_size,
    d = cond_size,
    tol = tol)
  
  return(prob)
}


