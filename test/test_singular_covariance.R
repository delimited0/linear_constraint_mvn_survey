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


d = 4
mu = muf(d)
Sigma = Sigmaf(d)
lb = lbf(d)
ub = ubf(d)
A = Af(d)

mu_t = as.numeric(A %*% mu)
lb_t = lb - mu_t
ub_t = ub - mu_t
Sigma_t = A %*% Sigma %*% t(A)


Q = chol(Sigma_t, pivot=TRUE)
pivot = attr(Q, "pivot")
t(Q[, order(pivot)]) %*% Q[, order(pivot)]
crossprod(Q[, order(pivot)])


true_val = (pnorm(2) - pnorm(-2))^d / 2
true_val

prob <- nvmix::pnvmix(ub_t, lb_t, qmix = "constant", 
                      loc = mu_t, factor = A, 
                      standardized = TRUE,
                      verbose = TRUE)
prob

nvmix::pnvmix(ub_t, lb_t, qmix = "constant", loc = mu_t, scale = Sigma_t,
              standardized = FALSE, control = list(cholesky.tol = 1))

nvmix:::cholesky_(Sigma_t, 1e-1)

mvtnorm::pmvnorm(lb_t, ub_t, mu_t, sigma = Sigma_t)

epmgpr::pmvn2(mu, Sigma, lb, ub, A)

epmgpr::pmvn(lb_t, ub_t, mu_t, Sigma_t)

met::pmvn(mu_t, Sigma_t, lb_t, ub_t)
  
foo = eigen(Sigma_t)

foo$values
