'fernandez_small

Usage:
  fernandez_small.R (--method_conf=<method_conf>) (--dim_conf=<dim_conf>) (--result_path=<result_path>) [--n_reps=n_reps] [--seed=seed] [--n_cores=n_cores] [--n_blas_threads=n_blas_threads]
  fernandez_small.R (-h|--help)

Options:
  -h --help  Usage.
  --method_conf=<method_conf>  Method configuration.
  --dim_conf=<dim_conf>  Dimension configuration.
  --result_path=<result_path>  Configuration specific sampler output.
  --n_reps=<n_reps>  Number of repetitions.
  --seed=seed  Seed.
  --n_cores=n_cores  Number of cores.
  --n_blas_threads=n_blas_threads  Number of BLAS threads.
' -> doc

opts = docopt::docopt(doc, version = 'fernandez_small 1.0')

method_conf = opts$method_conf
dim_conf = opts$dim_conf
result_path = opts$result_path
seed = as.numeric(opts$seed)
n_cores = as.numeric(opts$n_cores)
n_blas_threads = as.numeric(opts$n_blas_threads)
n_reps = as.integer(opts$n_reps)

# default arguments ------------------------------------------------------
if (is.null(seed)) 
  seed = 2022
if (is.null(n_cores))
  n_cores = 4
if (is.null(n_blas_threads))
  n_blas_threads = 1
if (is.null(n_reps)) 
  n_reps = 4

# hard coded arguments for debugging --------------------------------------
# method_conf = "experiments/fernandez_small/method_conf.json"
# dim_conf = "experiments/fernandez_small/dim_conf.json"
# result_path = "experiments/fernandez_small/test_results/"

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
RhpcBLASctl::blas_set_num_threads(n_blas_threads)
RhpcBLASctl::omp_set_num_threads(n_blas_threads)
doFuture::registerDoFuture()
future::plan(future::multicore, workers = n_cores)

# each job needs enough memory, 6GB
options(future.globals.maxSize = 6000 * 1024^2)

# start up summary --------------------------------------------------------
print(paste0("Running Fernandez 2007 small probability estimation with ",
             future::nbrOfWorkers(), " workers, each with",
             RhpcBLASctl::blas_get_num_procs(), " threads"))
print(paste0("Comparing ", nrow(methods), " methods:"))
sapply(methods, function(x) x$method)
print(paste0("Evaluating dimensions ", 
             paste0(dimensions, collapse = ", ")))
print(paste0(n_reps, " repetitions"))

# run simulation ----------------------------------------------------------
settings = expand.grid(method = methods, dimension = dimensions)
n_settings = nrow(settings)

progressr::handlers("progress")
progressr::with_progress({
  p = progressr::progressor(along = 1:(n_settings))
  
  for (method in methods) {
    
    # handle method output directories
    method_result_path = 
      paste0(result_path, "/", method$method, "/")
    if (!dir.exists(method_result_path)) 
      dir.create(method_result_path, recursive=TRUE)
    
    for (d in dimensions) {
      p(message = sprintf("Computing %s, dimension %d", method, d))
      
      # problem-dimension specific settings
      problem_params = list(
        mu = rep(0, d),
        Sigma = solve(.5 * diag(d) + .5 * rep(1, d) %*% t(rep(1, d))),
        lb = rep(.5, d),
        ub = rep(1, d)
      )
      
      # all input arguments
      params = c(problem_params, method$parameters)
      
      results = 
        foreach(i = 1:n_reps, .inorder = FALSE, .options.RNG = seed,
                .export = ls(globalenv()),
                .errorhandling = "remove") %dorng% #%dorng% 
        {
          tictoc::tic()
          result = do.call(method$method, params)
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


# progressr::handlers("progress")
# progressr::with_progress({
#   p = progressr::progressor(along = 1:(n_settings))
#   
#   foreach(i = 1:n_settings, .inorder = FALSE, .options.RNG = seed,
#           .export = ls(globalenv())) %dorng% 
#     {
#       method = settings[[i, "method"]]
#       d = settings[[i, "dimension"]]
#       
#       p(message = sprintf("%s, dimension %d", method, d))
#       
#       # problem-dimension specific settings
#       problem_params = list(
#         mu = rep(0, d),
#         Sigma = solve(.5 * diag(d) + .5 * rep(1, d) %*% t(rep(1, d))),
#         lb = rep(.5, d),
#         ub = rep(1, d)
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
#         paste0(result_path, "/", method$method, "/")
#       if (!dir.exists(method_result_path)) 
#         dir.create(method_result_path, recursive=TRUE)
#       
#       # save result
#       saveRDS(result, paste0(method_result_path, "d=", d))
#     }
# }, enable = TRUE)
