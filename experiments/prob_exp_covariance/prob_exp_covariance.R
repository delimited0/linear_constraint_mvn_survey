'prob_exp_covariance

Usage:
  prob_exp_covariance.R (--method_conf=<method_conf>) (--dim_conf=<dim_conf>) (--dep=<dep>) (--result_path=<result_path>) [--n_reps=n_reps] [--seed=seed] [--n_cores=n_cores] [--n_blas_threads=n_blas_threads]
  prob_exp_covariance.R (-h|--help)

Options:
  -h --help  Usage.
  --method_conf=<method_conf>  Method configuration.
  --dim_conf=<dim_conf>  Dimension configuration.
  --dep=<dep> Exponential covariance kernel dependence parameter.
  --result_path=<result_path>  Configuration specific sampler output.
  --n_reps=<n_reps>  Number of repetitions.
  --seed=seed  Seed.
  --n_cores=n_cores  Number of cores.
  --n_blas_threads=n_blas_threads  Number of threads for BLAS
' -> doc

opts = docopt::docopt(doc, version = 'prob_exp_covariance 1.0')

# Experiment from section 6.2, Genton et al 2018

method_conf = opts$method_conf
dim_conf = opts$dim_conf
dep = as.numeric(opts$dep)
result_path = opts$result_path
n_reps = opts$n_reps
seed = as.numeric(opts$seed)
n_cores = as.numeric(opts$n_cores)
n_blas_threads = as.numeric(opts$n_blas_threads)

# default arguments ------------------------------------------------------
if (is.null(seed)) 
  seed = 2022
if (is.null(n_cores)) 
  n_cores = 4
if (is.null(n_reps)) 
  n_reps = 4
if (is.null(n_blas_threads)) 
  n_blas_threads = 1

# hard coded arguments for debugging --------------------------------------
# method_conf = "experiments/prob_exp_covariance/method_conf.json"
# dim_conf = "experiments/prob_exp_covariance/test_dim_conf.json"
# result_path = "experiments/prob_exp_covariance/test_results/"
# dep = .1

# libraries ---------------------------------------------------------------
library(here)
library(doFuture)
library(progressr)
library(doRNG)
library(foreach)
library(tictoc)
library(fields)

source(here("prob_wrapper.R"))

# read settings -----------------------------------------------------------
methods = jsonlite::read_json(method_conf, simplifyVector = FALSE)
dimensions = jsonlite::read_json(dim_conf, simplifyVector = TRUE)

n_methods = length(methods)
n_dims = length(dimensions)

# setup output directories -------------------------------------------
if (!dir.exists(result_path)) dir.create(result_path)

# parallel set up ---------------------------------------------------------
RhpcBLASctl::blas_set_num_threads(n_blas_threads)  
RhpcBLASctl::omp_set_num_threads(n_blas_threads)
doFuture::registerDoFuture()
# future::plan(future::multisession, workers = n_cores)
future::plan(future::multicore, workers = n_cores)

# each job needs enough memory, 6GB
options(future.globals.maxSize = 6000 * 1024^2)

# start up summary --------------------------------------------------------
print(paste0("Running exponential covariance probability estimation with ",
             future::nbrOfWorkers(), " workers, each with",
             RhpcBLASctl::blas_get_num_procs(), " threads"))
print(paste0("Comparing ", nrow(methods), " methods:"))
sapply(methods, function(x) x$method)
print(paste0("Evaluating dimensions ", 
             paste0(dimensions, collapse = ", ")))
print(paste0(n_reps, " repetitions"))

# run simulation ----------------------------------------------------------

progressr::handlers("progress")
progressr::with_progress({
  p = progressr::progressor(along = 1:(n_dims*n_methods))
  
  for (method in methods) {
  
    # handle output directories
    method_result_path = 
      paste0(result_path, "/dep=", dep, "/", method$method, "/")
    if (!dir.exists(method_result_path)) 
      dir.create(method_result_path, recursive=TRUE)
    
    for (d in dimensions) {
      
      p(message = sprintf("Computing %s, dimension %d", method, d))
      
      set.seed(seed)
      locs = matrix(runif(2*d), ncol = 2)
      idx = tlrmvnmvt::zorder(locs)
      Sigma = exp(-rdist(locs) / dep)
      
      problem_params = list(
        mu = rep(0, d),
        Sigma = Sigma[idx, idx],
        lb = rep(-Inf, d),
        ub = rnorm(d, mean = 4, sd = .5)
      )
      
      # all input arguments
      params = c(problem_params, method$parameters)
      
      results = 
        foreach(i = 1:n_reps, .inorder = FALSE, .options.RNG = seed,
                .export = ls(globalenv()),
                .errorhandling = "remove") %dorng% #%dorng% 
        {
          # problem-dimension specific settings
          
          
          tictoc::tic()
          result = do.call(method$method, args = params)
          elapsed = tictoc::toc(quiet = TRUE)
          
          attr(result, "method") = method$method
          attr(result, "runtime") = elapsed$toc - elapsed$tic
          attr(result, "d") = d
          attr(result, "rep") = i
          
          return(result)
        }
      
      # save results together
      saveRDS(results, paste0(method_result_path, "d=", d))
    }
  }
}, enable = TRUE)


# settings = expand.grid(method = methods, dimension = dimensions)
# n_settings = nrow(settings)
# progressr::handlers("progress")
# progressr::with_progress({
#   p = progressr::progressor(along = 1:(n_settings))
#   
#   foreach(i = 1:n_settings, .inorder = FALSE, .options.RNG = seed,
#           .export = ls(globalenv()), .errorhandling = "remove") %dorng% 
#     {
#       method = settings[[i, "method"]]
#       d = settings[[i, "dimension"]]
#       
#       p(message = sprintf("%s, dimension %d", method, d))
#       
#       # problem-dimension specific settings
#       locs = matrix(runif(2*d), ncol = 2)
#       idx = tlrmvnmvt::zorder(locs)
#       Sigma = exp(-rdist(locs))
#       
#       problem_params = list(
#         mu = rep(0, d),
#         Sigma = Sigma[idx, idx],
#         lb = rep(-Inf, d),
#         ub = rnorm(d, mean = 4, sd = .5)
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
#         paste0(result_path, "/dep=", dep, "/", method$method, "/")
#       if (!dir.exists(method_result_path)) 
#         dir.create(method_result_path, recursive=TRUE)
#       
#       # save result
#       saveRDS(result, paste0(method_result_path, "d=", d))
#     }
# }, enable = TRUE)
