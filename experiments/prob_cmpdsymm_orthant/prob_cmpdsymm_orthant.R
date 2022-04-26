'prob_cmpdsymm_orthant

Usage:
  prob_cmpdsymm_orthant.R (--method_conf=<method_conf>) (--dim_conf=<dim_conf>) (--corr=<corr>) (--result_path=<result_path>) [--seed=seed] [--n_cores=n_cores]
  prob_cmpdsymm_orthant.R (-h|--help)

Options:
  -h --help  Usage.
  --method_conf=<method_conf>  Method configuration.
  --dim_conf=<dim_conf>  Dimension configuration.
  --corr=<corr> Correlation.
  --result_path=<result_path>  Configuration specific sampler output.
  --reps=<reps>  Number of repetitions.
  --seed=seed  Seed.
  --n_cores=n_cores  Number of cores.
' -> doc

opts = docopt::docopt(doc, version = 'prob_cmpdsymm_orthant 1.0')

method_conf = opts$method_conf
dim_conf = opts$dim_conf
corr = as.numeric(opts$corr)
result_path = opts$result_path
seed = as.numeric(opts$seed)
n_cores = as.numeric(opts$n_cores)

# default arguments ------------------------------------------------------
if (is.null(seed)) seed = 2022
if (is.null(n_cores)) n_cores = 1

# hard coded arguments for debugging --------------------------------------
# method_conf = "experiments/prob_cmpdsymm_orthant/method_conf.json"
# dim_conf = "experiments/prob_cmpdsymm_orthant/dim_conf.json"
# result_path = "experiments/prob_cmpdsymm_orthant/test_results/"

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

# setup output directories -------------------------------------------
if (!dir.exists(result_path)) dir.create(result_path)

# parallel set up ---------------------------------------------------------
RhpcBLASctl::blas_set_num_threads(1)  # no hyperthreading in BLAS
RhpcBLASctl::omp_set_num_threads(1)
doFuture::registerDoFuture()
future::plan(future::multicore, workers = n_cores)

# start up summary --------------------------------------------------------
print(paste0("Running compound symmetric orthant probability estimation with ",
             future::nbrOfWorkers(), " workers."))
print(paste0("Comparing ", nrow(methods), " methods:"))
sapply(methods, function(x) x$method)
print(paste0("Evaluating dimensions ", 
             paste0(dimensions, collapse = ", ")))

# run simulation ----------------------------------------------------------
settings = expand.grid(method = methods, dimension = dimensions)
n_settings = nrow(settings)

progressr::handlers("progress")
progressr::with_progress({
  p = progressr::progressor(along = 1:(n_settings))
  
  results = foreach(i = 1:n_settings, .inorder = FALSE, .options.RNG = seed,
          .export = ls(globalenv()), .errorhandling = "remove") %dorng% 
    {
      method = settings[[i, "method"]]
      d = settings[[i, "dimension"]]
      
      p(message = sprintf("%s, dimension %d", method, d))
      
      # problem-dimension specific settings
      problem_params = list(
        mu = rep(0, d),
        Sigma = corr * diag(d) + (1-corr) * rep(1, d) %*% t(rep(1, d)),
        lb = rep(0, d),
        ub = rep(Inf, d)
      )
      
      # method specific settings
      param_string = paste(
        names(method$parameters),
        method$parameters, 
        sep = "=", collapse = ", "
      )
      
      # all input arguments
      params = c(problem_params, method$parameters)
      
      tictoc::tic()
      result = do.call(method$method, params)
      elapsed = tictoc::toc(quiet = TRUE)
      
      attr(result, "method") = method$method
      attr(result, "runtime") = elapsed$toc - elapsed$tic
      
      # handle output directories
      method_result_path = 
        paste0(result_path, "/corr=", corr, "/", method$method, "/")
      if (!dir.exists(method_result_path)) 
        dir.create(method_result_path, recursive=TRUE)
      
      # save result
      saveRDS(result, paste0(method_result_path, "d=", d))
    }
}, enable = TRUE)


