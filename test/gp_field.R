source("prob_wrapper.R")

d = 512

locs = matrix(runif(2*d), ncol = 2)

idx = tlrmvnmvt::zorder(locs)

mu = rep(0, d)
Sigma = exp(-rdist(locs))
Sigma = Sigma[idx, idx]
lb = rep(-Inf, d)
ub = rnorm(d, mean = 4, sd = .5)


tlr(mu, Sigma, lb, ub, n_batch_mc = 1000)
sov(mu, Sigma, lb, ub, n_batch_mc = 1000)

met(mu, Sigma, lb, ub, n_batch_mc = 1000)
