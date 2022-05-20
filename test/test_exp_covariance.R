source("prob_wrapper.R")

library(fields)

d = 16
locs = matrix(runif(2*d), ncol = 2)
idx = tlrmvnmvt::zorder(locs)
Sigma = exp(-rdist(locs))

mu = rep(0, d)
Sigma = Sigma[idx, idx]
lb = rep(-Inf, d)
ub = rnorm(d, mean = 4, sd = .5)


hmvn_est = hmvn(mu, Sigma, lb, ub, n_batch_mc = 1000, block_size = 4, n_est = 10)
uvcdn_est = uvcdn(mu, Sigma, lb, ub)
hblkcdn_est = hblkcdn(mu, Sigma, lb, ub, cond_size = 4, block_size = 4)
lcg_est = lcg(mu, Sigma, lb, ub, 
              n_sub_samples = 16,
              n_sub_skip = 10,
              n_hdr_samples = 420,
              n_hdr_skip = 2)

met_est = met(mu, Sigma, lb, ub, n_batch_mc = 1000)


problem_params = list(
  mu = mu, 
  Sigma = Sigma,
  lb = lb,
  ub = ub
)

hmvn_params = c(problem_params, list(block_size = 4, n_batch_mc = 1000))
do.call("hmvn", hmvn_params)
