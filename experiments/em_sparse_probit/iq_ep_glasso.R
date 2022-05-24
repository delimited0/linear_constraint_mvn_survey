'iq_ep_glasso

Usage:
  iq_ep_glasso.R (--output_path=<output_path>) (--n_reps=<n_reps>) [--max_iter=max_iter] [--seed=seed] [--n_cores=n_cores] [--n_blas_threads=n_blas_threads]
  iq_glasso.R (-h|--help)

Options:
  -h --help  Usage.
  --output_path=<output_path>  Configuration specific sampler output.
  --n_reps=<reps>  Number of repetitions.
  --max_iter=max_iter  Maximum EM iterations.
  --seed=seed  Seed.
  --n_cores=n_cores  Number of cores.
  --n_blas_threads=n_blas_threads  Number of threads for multithread BLAS
' -> doc

opts = docopt::docopt(doc, version = 'EP sparse probit 1.0')

output_path = opts$output_path
seed = as.numeric(opts$seed)
n_cores = as.numeric(opts$n_cores)
n_blas_threads = as.numeric(opts$n_blas_threads)
n_reps = as.numeric(opts$n_reps)
max_iter = as.numeric(opts$max_iter)

# default arguments -------------------------------------------------------
if (is.null(seed)) 
  seed = 2022
if (is.null(n_cores)) 
  n_cores = 1
if (is.null(n_blas_threads)) 
  n_blas_threads = 1
if (is.null(max_iter)) 
  max_iter = 200

# hard coded arguments for debugging --------------------------------------
# output_path = "experiments/em_sparse_probit/test_output"
# problem_idx = 1
# j = 1
# max_iter = 300
# n_reps = 4

# libraries ---------------------------------------------------------------
library(here)
library(doFuture)
library(progressr)
library(doRNG)
library(foreach)
library(tictoc)

source(here("experiments", "em_sparse_probit", "em_glasso_probit.R"))

# setup output directories -------------------------------------------
if (!dir.exists(output_path)) dir.create(output_path)

# parallel set up ---------------------------------------------------------
RhpcBLASctl::blas_set_num_threads(n_blas_threads)
RhpcBLASctl::omp_set_num_threads(n_blas_threads)
doFuture::registerDoFuture()
# future::plan(future::multisession, workers = n_cores)
future::plan(future::multicore, workers = n_cores)

# start up summary --------------------------------------------------------
print(paste0("Running EP-EM sparse covariance probit estimation 
             small examples with ",
             future::nbrOfWorkers(), " workers, ", 
             n_blas_threads, " BLAS threads / worker."))

print(sprintf("%d repetitions", n_reps))

handlers("progress")

# run simulation ----------------------------------------------------------

n_obs = 2000
d = 1
n_mc_samples = 50

true_precisions = list(
  "3d" = matrix(c(1, .5, 0, .5, 1, 0, 0, 0, 1), nrow = 3),
  "5d_a" = matrix(c(1, 0,    0,    0, 0, 
                    0, 1,    0,    0, 0,
                    0, 0,    1, -.53, 0,
                    0, 0, -.53,    1, 0,
                    0, 0,    0,    0, 1), nrow = 5),
  "5d_b" = matrix(c(  1, .06,   0, .37,   0,
                      .06,   1,   0, .16,   0,
                      0,   0,   1,   0, .56,
                      .37, .16,   0,   1,   0,
                      0,   0, .56,   0,   1), nrow = 5)
)

for (problem_idx in 1:length(true_precisions)) {
  
  experiment_name = names(true_precisions)[problem_idx]
  print(paste0("Experiment: ", experiment_name))
  
  # true parameter values
  coef_true = rep(0, d)
  Prec_true = true_precisions[[problem_idx]]
  Sigma_true = solve(Prec_true)
  n_choices = nrow(Sigma_true)
  
  Sigma_init = diag(n_choices)
  beta_init = rep(0, d)
  
  # regularization amount
  penalty = sqrt(log(n_choices) / n_obs) * .1
  
  method_output_path = paste0(output_path, "/", experiment_name, "/EP/")
  if (!dir.exists(method_output_path)) 
    dir.create(method_output_path, recursive = TRUE)
  
  progressr::with_progress({
    progger = progressr::progressor(along = 1:n_reps)    
    
    foreach(i = 1:n_reps, .inorder = FALSE, .options.RNG = seed,
            .errorhandling = "remove") %dorng% {
              
              # simulate data
              X = lapply(1:n_obs, function(obs) {
                matrix(runif(n_choices * d, min = -.5, max = .5),
                       nrow = n_choices, ncol = d)
              })
              Z = t(sapply(X, function(x) {
                mvtnorm::rmvnorm(1, x %*% coef_true, Sigma_true) 
              }))
              Y = t(apply(Z, 1, function(x) x >= x[which.max(x)])) 
              
              tic()
              result = probit_ep_glasso(
                X = X, Y = Y, 
                penalty = penalty, 
                Sigma_init = Sigma_init, beta_init = beta_init, 
                corr = TRUE, pen_diag = FALSE, 
                max_iter = max_iter)
              elapsed = toc(quiet=TRUE)
              
              attr(result, "method") = "EP"
              attr(result, "runtime") = elapsed$toc - elapsed$tic
              attr(result, "rep") = i
              
              saveRDS(result, paste0(method_output_path, "rep=", i))
              
              progger(sprintf("Problem %s, method: %s, rep=%i", 
                              experiment_name, "EP", i))
            }
  }, enable = TRUE)
}


