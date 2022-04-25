source("prob_wrapper.R")

n = 1000
d = 16
mu = rep(0, d)
Sigma = .5*diag(d) + .5*rep(1, d) %*% t(rep(1, d))
lb = rep(0, d)
ub = rep(Inf, d)

sov(mu, Sigma, lb, ub, n_batch_mc = n)

tlr(mu, Sigma, lb, ub, n_batch_mc = n)

epmgp(mu, Sigma, lb, ub)

ghk(mu, Sigma, lb, ub, n_batch_mc = n)

lcg(mu, Sigma, lb, ub, 
    n_sub_samples = 16, n_sub_skip = 10,
    n_hdr_samples = 420, n_hdr_skip = 2)

met(mu, Sigma, lb, ub, n_batch_mc = n, n_est = 10)

bvcdn(mu, Sigma, lb, ub)


# Fernandez 2007 ----------------------------------------------------------
d = 40
mu = rep(0, d)
Sigma = solve(.5 * diag(d) + .5 * rep(1, d) %*% t(rep(1, d)))
lb = rep(.5, d)
ub = rep(1, d)
n = 1000

tlr(mu, Sigma, lb, ub, n_batch_mc = n)

met(mu, Sigma, lb, ub, n_batch_mc = n, n_est = 10)
TruncatedNormal::pmvnorm(mu, Sigma, lb, ub, tyep = "qmc")


# special conditioning tests ----------------------------------------------

d = 64
m = 64

mu = rep(0, d)
Sigma = .5*diag(d) + .5*rep(1, d) %*% t(rep(1, d))
lb = rep(0, d)
ub = rep(Inf, d)

bvcdn(mu, Sigma, lb, ub)

hblkcdn(mu, Sigma, lb, ub, cond_size = 2, block_size = )#  patrick stinks