n = 1000
d = 2
mu = rep(0, d)
Sigma = .5*diag(d) + .5*rep(1, d) %*% t(rep(1, d))
lb = rep(0, d)
ub = rep(Inf, d)


gibbs_ry2004(n, mu, Sigma, lb, ub, initial = c(1, 1), burnin = 100, tuvn_sampler = "lg2015")

epess_samples = epess(n, mu, Sigma, lb, ub, initial = c(1, 1), N=1, J=1, burnin = 100)
plot(epess_samples)

sampletmvn::sample_gibbs_ry2004(n, mu, Sigma, lb, ub, A = NULL,init = c(1, 1), tuvn_sampler = "lg2015", burn = 100, thin = 0)
