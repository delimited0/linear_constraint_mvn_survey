library(here)

source(here("sampling_wrapper.R"))
source(here("experiments", "em_sparse_probit", "em_glasso_probit.R"))

n_obs = 2000
p = 1
n_mc_samples = 50

Prec_true = matrix(c(1, .5, 0, .5, 1, 0, 0, 0, 1), nrow = 3)

# true parameter values
coef_true = rep(0, p)
Prec_true = true_precisions[[problem_idx]]
Sigma_true = solve(Prec_true)
n_choices = nrow(Sigma_true)

Sigma_init = diag(n_choices)
beta_init = rep(0, p)

# regularization amount
penalty = sqrt(log(n_choices) / n_obs) * .1

# simulate data
X = lapply(1:n_obs, function(i) {
  matrix(runif(n_choices * p, min = -.5, max = .5),
         nrow = n_choices, ncol = p)
})
Z = t(sapply(X, function(x) {
  mvtnorm::rmvnorm(1, x %*% coef_true, Sigma_true) 
}))
Y = t(apply(Z, 1, function(x) x >= x[which.max(x)])) 

lb = rep(0, n_choices-1)
ub = rep(Inf, n_choices-1)

result = probit_mcmc_glasso(
  X = X, Y = Y, 
  sampler = method, sampler_params = params,
  penalty = penalty, 
  Sigma_init = Sigma_init, beta_init = beta_init, 
  corr = FALSE, pen_diag = FALSE, 
  n_mc_samples = n_mc_samples)

i = 1
A = max_polytope(Y[i, ])
initial = 2*Y[i, ] - 1
sampletmvn::sample_gibbs_lg2015(n = n_mc_samples, 
                    mu = X[[i]] %*% beta_init,
                    Sigma = Sigma_init, 
                    lower = lb, upper = ub, A = A, 
                    init = initial,
                    burn = 0,
                    thin = 1,
                    tuvn_sampler = "lg2015")  
