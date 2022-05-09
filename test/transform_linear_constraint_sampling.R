muf <- function(d) rep(0, d)
Sigmaf <- function(d) diag(d)
lbf <- function(d) rep(-Inf, 2*d)
ubf <- function(d) c(0, rep(2, 2*d-1))
Af <- function(d) {
  lower_bounds <- -diag(d)
  upper_bounds <- diag(d)
  upper_bounds[1, ] <- c(2, 1, rep(0, d-2))
  A <- rbind(upper_bounds, lower_bounds)
  return(A)
}

d = 2
mu = muf(d)
Sigma = Sigmaf(d)
lb = lbf(d)
ub = ubf(d)
A = Af(d)

init = rep(-1, d)

source("sampling_wrapper.R")

# direct polytope sampling ------------------------------------------------
n = 1000

direct_sampes = gibbs_lg2015(n = n, mu = mu, Sigma = Sigma, lb = lb, ub = ub, 
             A = A, initial = init, tuvn_sampler = "lg2015")


# MP inverse transform sampling -------------------------------------------

G = MASS::ginv(A)

G %*% A  # G is left inverse for A

lb_t = lb - A %*% mu
ub_t = ub - A %*% mu
Sigma_t = A %*% Sigma %*% t(A)
mu_t = as.numeric(mu - A %*% mu)

# unless you can handle singular covariance this isn't going to work
trans_samples = gibbs_cov(n, mu_t, Sigma_t, lb_t, ub_t, initial = init)


