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
d_t = sqrt( 1 / (diag(Sigma_t)) )

lb_f = d_t * lb_t
ub_f = d_t * ub_t
Corr_t = cov2cor(Sigma_t)


# why is cholesky failing on Sigma_t? ----
cholesky_ <- function(mat, tol = 1e-12)
{
  n <- dim(mat)[1] # dimension
  stopifnot(dim(mat)[2] == n)
  ## First try 'chol()' (=>fast)
  C <- tryCatch(t(chol(mat)), error = function(e) e)
  if(is.matrix(C) && all.equal(dim(C), rep(n, 2))) {
    ## C is the desired Cholesky factor
    ## Grab diagonal
    diag.elements <- diag(C)
  } else {
    ## In this case, 't(chol(scale))' returned an error so that we manually
    ## compute the Cholesky factor of the *singular* matrix 'mat'
    C <- matrix(0, ncol = n, nrow = n) # initialize Cholesky factor
    diag.elements <- rep(NA, n)
    for(col in 1:n) {
      dsq <- mat[col, col] - sum(C[col, 1:(col-1)]^2) # C(col,col)^2
      browser()
      if(dsq < 0) stop("Matrix not positive semi-definite")
      d <- if(dsq < abs(mat[col, col] * tol)) 0 else sqrt(dsq) # set 0 if too small
      C[col, col] <- d
      diag.elements[col] <- d # store diagnonal element
      if(col < n && d > 0) { # calculate the remaining elements in column 'col'
        for(row in (col+1):n) {
          C[row, col] <- (mat[row, col] -
                            sum(C[row, 1:(row-1)]*C[col, 1:(row-1)]))/d
        }
      }
    }
  }
  list(C = C, D = diag.elements)
}

cholesky_(Sigma_t)

Q = chol(Sigma_t, pivot=TRUE)
pivot = attr(Q, "pivot")
t(Q[, order(pivot)]) %*% Q[, order(pivot)]
crossprod(Q[, order(pivot)])


true_val = (pnorm(2) - pnorm(-2))^d / 2
true_val

prob <- nvmix::pnvmix(ub_t, lb_t, qmix = "constant", 
                      loc = mu_t, factor = A, 
                      standardized = TRUE,
                      verbose = TRUE, 
                      control = list(B = 10, max.iter.rqmc = 10, 
                                     increment = "num.init",
                                     fun.eval = c(1000, 1000)))
prob

nvmix::pnvmix(ub_t, lb_t, qmix = "constant", loc = mu_t, scale = Corr_t,
              standardized = TRUE, control = list(cholesky.tol = 1e-1))

nvmix:::cholesky_(Sigma_t, 1e-1)
nvmix:::cholesky_(Corr_t, 1e-1)

eigen(Corr_t)$values

mvtnorm::pmvnorm(lb_t, ub_t, mu_t, sigma = Sigma_t)

epmgpr::pmvn2(mu, Sigma, lb, ub, A)

epmgpr::pmvn(lb_t, ub_t, mu_t, Sigma_t)

met::pmvn(mu_t, Sigma_t, lb_t, ub_t)
  
foo = eigen(Sigma_t)

foo$values



