'prob_cmpdsymm_orthant

Usage:
  prob_cmpdsymm_orthant.R (--method_conf=<method_conf>) (--dim_conf=<dim_conf>) (--corr=<corr>) (--result_path=<result_path>) [--n_reps=n_reps] [--seed=seed] [--n_cores=n_cores] [--n_blas_threads=n_blas_threads]
  prob_cmpdsymm_orthant.R (-h|--help)

Options:
  -h --help  Usage.
  --method_conf=<method_conf>  Method configuration.
  --dim_conf=<dim_conf>  Dimension configuration.
  --corr=<corr> Correlation.
  --result_path=<result_path>  Configuration specific sampler output.
  --n_reps=<reps>  Number of repetitions.
  --seed=seed  Seed.
  --n_cores=n_cores  Number of cores.
  --n_blas_threads=n_blas_threads  Number of BLAS threads.
' -> doc

opts = docopt::docopt(doc, version = 'prob_cmpdsymm_orthant 1.0')

method_conf = opts$method_conf
dim_conf = opts$dim_conf
corr = as.numeric(opts$corr)
result_path = opts$result_path
n_reps = opts$n_reps
seed = as.numeric(opts$seed)
n_cores = as.numeric(opts$n_cores)
n_blas_threads = as.numeric(opts$n_blas_threads)

# default arguments ------------------------------------------------------
if (is.null(seed)) seed = 2022
if (is.null(n_cores)) n_cores = 1
if (is.null(n_blas_threads)) n_blas_threads = 1

# hard coded arguments for debugging --------------------------------------
# method_conf = "experiments/prob_cmpdsymm_orthant/method_conf.json"
# dim_conf = "experiments/prob_cmpdsymm_orthant/test_dim_conf.json"
# result_path = "experiments/prob_cmpdsymm_orthant/test_results/"
# corr = .5

# libraries ---------------------------------------------------------------
library(here)
library(doFuture)
library(progressr)
library(doRNG)
library(foreach)
library(tictoc)

source(here("prob_wrapper.R"))

# read settings -----------------------------------------------------------
methods = jsonlite::read_json(method_conf, simplifyVector = FALSE)
dimensions = jsonlite::read_json(dim_conf, simplifyVector = TRUE)

n_methods = length(methods)
n_dims = length(dimensions)

# setup output directories -------------------------------------------
if (!dir.exists(result_path)) dir.create(result_path)

# parallel set up ---------------------------------------------------------
RhpcBLASctl::blas_set_num_threads(n_blas_threads)  # no hyperthreading in BLAS
RhpcBLASctl::omp_set_num_threads(n_blas_threads)
doFuture::registerDoFuture()
# future::plan(future::multicore, workers = n_cores)
future::plan(future::multisession, workers = n_cores)

# each job needs enough memory, 6GB
options(future.globals.maxSize = 10000 * 1024^2)

# start up summary --------------------------------------------------------
print(paste0("Running compound symmetric orthant probability estimation with ",
             future::nbrOfWorkers(), " workers, each with ",
             RhpcBLASctl::blas_get_num_procs(), " threads"))
print(paste0("Comparing ", n_methods, " methods:"))
sapply(methods, function(x) x$method)
print(paste0("Evaluating dimensions ", 
             paste0(dimensions, collapse = ", ")))
print(paste0(n_reps, " repetitions"))

# run simulation ----------------------------------------------------------
# settings = expand.grid(method = methods, dimension = dimensions)
# n_settings = nrow(settings)

handlers("progress")
with_progress({
  p = progressor(along = 1:(n_dims*n_methods))
  
  for (method in methods) {
    
    # handle output directories
    method_result_path = 
      paste0(result_path, "/corr=", corr, "/", method$method, "/")
    if (!dir.exists(method_result_path)) 
      dir.create(method_result_path, recursive=TRUE)
    
    for (d in dimensions) {
      
      p(message = sprintf("Computing %s, dimension %d", method, d))
      
      # problem-dimension specific settings
      problem_params = list(
        mu = rep(0, d),
        Sigma = corr * diag(d) + (1-corr) * rep(1, d) %*% t(rep(1, d)),
        lb = rep(0, d),
        ub = rep(Inf, d)
      )
      
      # all input arguments
      params = c(problem_params, method$parameters)
      
      # r = progressor(along = 1:n_reps)
      results = 
        foreach(i = 1:n_reps, .inorder = FALSE, .options.RNG = seed,
                .export = ls(globalenv()),
                .errorhandling = "remove") %dorng% 
        {
          # r(message = sprintf("Repetition %d", d))
          
          tictoc::tic()
          result = do.call(method$method, params)
          elapsed = tictoc::toc(quiet = TRUE)
          
          attr(result, "method") = method$method
          attr(result, "runtime") = elapsed$toc - elapsed$tic
          attr(result, "d") = d
          attr(result, "rep") = i
          
          # save result by repetition
          # saveRDS(results, paste0(method_result_path, "d=", d, "_rep=", i))
          return(result)
        }
      
      # save results together
      saveRDS(results, paste0(method_result_path, "d=", d))
    }  
  }
}, enable=TRUE)
  
  
#   results = foreach(i = 1:n_settings, .inorder = FALSE, .options.RNG = seed,
#           .export = ls(globalenv()), .errorhandling = "remove") %dorng% 
#     {
#       method = settings[[i, "method"]]
#       d = settings[[i, "dimension"]]
#       
#       p(message = sprintf("%s, dimension %d", method, d))
#       
#       # problem-dimension specific settings
#       problem_params = list(
#         mu = rep(0, d),
#         Sigma = corr * diag(d) + (1-corr) * rep(1, d) %*% t(rep(1, d)),
#         lb = rep(0, d),
#         ub = rep(Inf, d)
#       )
#       
#       # method specific settings
#       param_string = paste(
#         names(method$parameters),
#         method$parameters, 
#         sep = "=", collapse = ", "
#       )
#       
#       # all input arguments
#       params = c(problem_params, method$parameters)
#       
#       tictoc::tic()
#       result = do.call(method$method, params)
#       elapsed = tictoc::toc(quiet = TRUE)
#       
#       attr(result, "method") = method$method
#       attr(result, "runtime") = elapsed$toc - elapsed$tic
#       attr(result, "d") = d
#       
#       # handle output directories
#       method_result_path = 
#         paste0(result_path, "/corr=", corr, "/", method$method, "/")
#       if (!dir.exists(method_result_path)) 
#         dir.create(method_result_path, recursive=TRUE)
#       
#       # save result
#       saveRDS(result, paste0(method_result_path, "d=", d))
#     }
# }, enable = TRUE)


