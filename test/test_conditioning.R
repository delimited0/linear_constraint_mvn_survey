d = 64
mu = rep(0, d)
Sigma = .5*diag(d) + .5*rep(1, d) %*% t(rep(1, d))
lb = rep(0, d)
ub = rep(Inf, d)
ub_approx = rep(10, d)

# Mendell Elston univariate conditioning ----
odr = hccmvn::recurUniandBlkcmb(Sigma, lb, ub, d)+1

mu_odr = mu[odr]
lb_odr = lb[odr]
ub_odr = ub[odr]
Sigma_odr = Sigma[odr, odr]

prob = hccmvn::hccmvn(covM = Sigma_odr, a = lb_odr, b = ub_approx, m = d, d = 1, tol = 1e-4)
prob

prob = hccmvn::hccmvn(covM = Sigma_odr, a = lb_odr, b = ub_approx, m = d, d = 1, tol = 1e-4)
prob

# bivariate conditioning
prob = hccmvn::hccmvn(covM = Sigma, a = lb, b = ub, m = d, d = 2, tol = 1e-4)
prob


prob = hccmvn::hccmvn(covM = Sigma, a = lb, b = ub_approx, m = d, d = 2, tol = 1e-4)
prob

prob = hccmvn::hccmvn(covM = Sigma_odr, a = lb_odr, b = ub_odr, m = d, d = 2, tol = 1e-4)
prob

# q variate conditioning
prob = hccmvn::hccmvn(covM = Sigma, a = lb, b = ub_approx, m = d, d = 6, tol = 1e-4)
prob

# hierarchical block conditioning
odr = hccmvn::recurUniandBlkcmb(Sigma, lb, ub, 4)

mu_odr = mu[odr]
lb_odr = lb[odr]
ub_odr = ub[odr]
Sigma_odr = Sigma[odr, odr]

prob = hccmvn::hccmvn(covM = Sigma_odr, a = lb_odr, b = ub_approx, m = 32, d = 4, tol = 1e-4)
prob

prob = hccmvn::hccmvn(covM = Sigma, a = lb, b = ub_approx, m = 32, d = 4, tol = 1e-4)
prob



# a spatial exponential covariance example --------------------------------
library(fields)
source("prob_wrapper.R")

d = 512

locs = matrix(runif(2*d), ncol = 2)
idx = tlrmvnmvt::zorder(locs)

mu = rep(0, d)
Sigma = exp(-rdist(locs))
Sigma = Sigma[idx, idx]
lb = rep(-Inf, d)
lb_approx = rep(-10, d)
ub = rnorm(d, mean = 4, sd = .5)

odr = hccmvn::recurUniandBlkcmb(Sigma, lb_approx, ub, bsz_ = 32)+1

# univariate conditioning
prob = hccmvn::hccmvn(covM = Sigma, a = lb_approx, b = ub, m = d, d = 1, tol = 1e-4)
prob

# bivariate conditioning
prob = hccmvn::hccmvn(covM = Sigma, a = lb_approx, b = ub, m = d, d = 2, tol = 1e-4)
prob

# q variate conditioning
prob = hccmvn::hccmvn(covM = Sigma, a = lb_approx, b = ub, m = d, d = 4, tol = 1e-4)
prob

tlr(mu, Sigma, lb, ub, n_batch_mc = 1000)


# non square sizes --------------------------------------------------------
source("prob_wrapper.R")
d = 20
mu = rep(0, d)
Sigma = .5*diag(d) + .5*rep(1, d) %*% t(rep(1, d))
lb = rep(0, d)
ub = rep(Inf, d)

uvcdn(mu, Sigma, lb, ub)
bvcdn(mu, Sigma, lb, ub)
dvcdn(mu, Sigma, lb, ub, cond_size = 4)


# univariate reordering safety --------------------------------------------

source("prob_wrapper.R")
d = 10
mu = rep(0, d)
Sigma = .5*diag(d) + .5*rep(1, d) %*% t(rep(1, d))
lb = rep(0, d)
ub = rep(Inf, d)

uvcdn(mu, Sigma, lb, ub)
bvcdn(mu, Sigma, lb, ub)


