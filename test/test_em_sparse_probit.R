library(here)

source(here("sampling_wrapper.R"))
source(here("experiments", "em_sparse_probit", "em_glasso_probit.R"))

n_obs = 2000
p = 1
n_mc_samples = 5

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
Y = t(
  apply(Z, 1, function(z) {
    choice = rep(FALSE, n_choices)
    choice[which.max(z)] = TRUE
    return(choice)
  })
) 

result_ep = probit_ep_glasso(
  X = X, Y = Y, 
  penalty = penalty, 
  Sigma_init = Sigma_init, beta_init = beta_init, 
  corr = TRUE, pen_diag = TRUE, 
  max_iter = max_iter)

result_ep$Prec
result_ep$Sigma

plot(result_ep$glasso_llik)

prec_history = sapply(result_ep$Prec_history, function(P) P[1,2])
plot(
  sapply(result_ep$Prec_history, function(P) P[1,2])
)

result = probit_mcmc_glasso(
  X = X, Y = Y, 
  sampler = method, sampler_params = params,
  penalty = penalty, 
  Sigma_init = Sigma_init, beta_init = beta_init, 
  corr = FALSE, pen_diag = FALSE, 
  n_mc_samples = n_mc_samples)


# simulation test
i = 1
lb = rep(0, n_choices-1)
ub = rep(Inf, n_choices-1)
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

G = MASS::ginv(A)

Gt = MASS::ginv(t(A))
G %*% A

# a larger precision entry ------------------------------------------------

n_obs = 2000
p = 1
n_mc_samples = 5

# Prec_true =  10 * matrix(c(1, 0,    0,    0, 0, 
#                      0, 1,    0,    0, 0,
#                      0, 0,    1,   -.99, 0,
#                      0, 0,   -.99,    1, 0,
#                      0, 0,    0,    0, 1), nrow = 5)
Prec_true = solve(.5 * diag(5) + .5 * rep(1, 5) %*% t(rep(1, 5)))

# true parameter values
coef_true = rep(0, p)
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
Y = t(
  apply(Z, 1, function(z) {
    choice = rep(FALSE, n_choices)
    choice[which.max(z)] = TRUE
    return(choice)
  })
) 

result_ep = probit_ep_glasso(
  X = X, Y = Y, 
  penalty = penalty, 
  Sigma_init = Sigma_init, beta_init = beta_init, 
  corr = FALSE, pen_diag = FALSE, 
  max_iter = max_iter)

result_ep$Prec
result_ep$Sigma

plot(result_ep$glasso_llik)

prec_history = sapply(result_ep$Prec_history, function(P) P[3,4])
plot(
  sapply(result_ep$Prec_history, function(P) P[3,4])
)
plot(result_ep$glasso_llik, prec_history, type = "l")

result_mcmc = probit_mcmc_glasso(
  X = X, Y = Y, 
  sampler = "rsm", sampler_params = list(),
  penalty = penalty, 
  Sigma_init = Sigma_init, beta_init = beta_init, 
  corr = TRUE, pen_diag = TRUE, 
  n_mc_samples = n_mc_samples, max_iter = max_iter
)

result_mcmc$Prec
plot(result_mcmc$glasso_llik)

prec_history = sapply(result_mcmc$Prec_history, function(P) P[3,4])
plot(prec_history)

plot(result_mcmc$glasso_llik, prec_history, type = "l")
