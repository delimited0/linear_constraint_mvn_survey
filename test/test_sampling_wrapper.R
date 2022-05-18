
# cmpdsymm positive orthant  ----------------------------------------------


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

ehmc_samples = ehmc(n, mu, Sigma, lb, ub, initial = c(1, 1))
plot(ehmc_samples)


# trapezoid ---------------------------------------------------------------



# multinomial probit sparse covariance ------------------------------------
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

y = c(0, 1, 0)
A = max_polytope(y)

lb = rep(0, 2)
ub = rep(Inf, 2)

n = 1000

mu = rep(0, 3)
Sigma = diag(3)

ehmc_samples = ehmc(n, mu, Sigma, lb, ub, A = A, initial =  2*y - 1)
par(mfrow = c(2, 2))
plot(ehmc_samples[, c(1, 2)])
plot(ehmc_samples[, c(2, 3)])
plot(ehmc_samples[, c(1, 3)])

met_samples = met_mp(n, mu, Sigma, lb, ub, A, initial = 2*y - 1)
par(mfrow = c(2, 2))
plot(met_samples[, c(1, 2)])
plot(met_samples[, c(2, 3)])
plot(met_samples[, c(1, 3)])

liness_samples = liness(n, mu, Sigma, lb, ub, A, initial = 2*y - 1)
par(mfrow = c(2, 2))
plot(liness_samples[, c(1, 2)])
plot(liness_samples[, c(2, 3)])
plot(liness_samples[, c(1, 3)])
