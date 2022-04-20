'sample_symmetric_orthant

Usage:
  sample_symmetric_orthant.R (--method_conf=<method_conf>) (--dim_conf=<dim_conf>)  (--sample_path=<sample_path>) [--seed=seed] [--n_threads=n_threads] 
  sample_symmetric_orthant.R (-h|--help)

Options:
  -h --help  Usage.
  --method_conf=<method_conf>  Method configuration.
  --dim_conf=<dim_conf>  Dimension configuration.
  --sample_path=<sample_path>  Configuration specific sampler output.
  --seed=seed  Seed.
  --n_threads=n_threads  Number of cores.
' -> doc

opts = docopt::docopt(doc, version = 'sample_symmetric_orthant 1.0')

method_conf = opts$method_conf
dim_conf = opts$dim_conf
sample_path = opts$sample_path
seed = as.numeric(opts$seed)
n_threads = as.numeric(opts$n_threads)

# default arguments -------------------------------------------------------
if (is.null(seed)) seed = 2022
if (is.null(n_threads)) n_threads = 1

# hard coded arguments for debugging --------------------------------------
# method_conf = "experiments/sample_symmetric_orthant/test_method_conf.json"
# dim_conf = "experiments/sample_symmetric_orthant/test_dim_conf.json"
# sample_path = "experiments/sample_symmetric_orthant/samples"
# timing_path = "experiments/sample_symmetric_orthant/timings"

# libraries ---------------------------------------------------------------
library(here)
library(doFuture)
library(progressr)
library(doRNG)
library(foreach)

source(here("sampling_wrapper.R"))

# read settings ---------------------------------------------------------
methods = jsonlite::read_json(method_conf, simplifyVector = FALSE)
dimensions = jsonlite::read_json(dim_conf, simplifyVector = TRUE)

# setup output directories -------------------------------------------
if (!dir.exists(sample_path)) dir.create(sample_path)

# parallel set up ---------------------------------------------------------
RhpcBLASctl::blas_set_num_threads(1)  # no hyperthreading in BLAS
doFuture::registerDoFuture()
future::plan(future::multisession, workers = n_threads)
print(paste0("Running sampling test 1 with ", future::nbrOfWorkers(), " workers."))

# start up summary --------------------------------------------------------
print(paste0("Running symmetric orthant sampling with ",
             future::nbrOfWorkers(), " workers."))
print(paste0("Comparing ", nrow(methods), " methods:"))
sapply(methods, function(x) x$method)
print(paste0("Evaluating dimensions ", 
             paste0(dimensions, collapse = ", ")))

# run simulation ----------------------------------------------------------
settings = expand.grid(method = methods, dimension = dimensions)
n_settings = nrow(settings)

progressr::with_progress({
  p = progressr::progressor(along = 1:(n_settings))
  
  foreach(i = 1:n_settings, .inorder = FALSE, .options.RNG = seed,
          .export = ls(globalenv())) %dorng% 
    {
      
      p(sprintf("i=%g", i))
      
      method = settings[[i, "method"]]
      d = settings[[i, "dimension"]]
      
      # problem-dimension specific settings
      problem_params = list(
        n = 5000,
        mu = rep(0, d),
        Sigma = .5*diag(d) + .5*rep(1, d) %*% t(rep(1, d)),
        lb = rep(0, d),
        ub = rep(Inf, d),
        initial = rep(1, d)
      )
      
      # method specific settings
      param_string = paste(
        names(method$parameters),
        method$parameters, 
        sep = "=", collapse = ", "
      )
      
      # print( paste0("--- ", method$method, ", Dimension: ", d, "---") )
      # print( paste("Parameters: ", param_string) )
      
      # all input arguments
      params = c(problem_params, method$parameters)
      
      tictoc::tic()
      samples = do.call(method$method, params)
      elapsed = tictoc::toc(quiet = TRUE)
      
      attr(samples, "method") = method$method
      attr(samples, "runtime") = elapsed$toc - elapsed$tic
      
      # handle output directories
      method_sample_path = 
        paste0(sample_path, "/", method$method, "/")
      if (!dir.exists(method_sample_path)) 
        dir.create(method_sample_path, recursive=TRUE)
      
      # save samples
      saveRDS(samples, paste0(method_sample_path, "d=", d))
    } 
})






