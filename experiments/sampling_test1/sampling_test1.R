'sampling_test1.

Usage:
  sampling_test1.R (--method_conf=<method_conf>) (--dim_conf=dim_conf) (--sample_path=<sample_path>) (--timing_path=<timing_path>) [--seed=seed] [--n_threads=n_threads]
  sampling_test1.R (-h|--help)

' -> doc

args = docopt::docopt(doc, version = 'sampling_test1 1.0')

method_conf = args$method_conf
dim_conf = args$dim_conf
sample_path = args$sample_path
timing_path = args$timing_path
seed = as.numeric(args$seed)
n_threads = as.numeric(args$n_threads)

# default arguments -------------------------------------------------------
if (is.null(seed)) seed = 2022
if (is.null(n_threads)) n_threads = 1


# libraries ---------------------------------------------------------------
source("sampling_wrapper.R")

library(doFuture)
library(progressr)
library(doRNG)
library(foreach)

# hard coded arguments for debugging --------------------------------------
# method_conf = "experiments/sampling_test1/method_conf.json"
# dim_conf = "experiments/sampling_test1/dim_conf.json"

# read settings ---------------------------------------------------------
methods = jsonlite::read_json(method_conf, simplifyVector = FALSE)
dimensions = jsonlite::read_json(dim_conf, simplifyVector = TRUE)

# setup output directories -------------------------------------------
if (!dir.exists(sample_path)) dir.create(sample_path)
if (!dir.exists(timing_path)) dir.create(timing_path)


# parallel set up ---------------------------------------------------------
RhpcBLASctl::blas_set_num_threads(1)  # no hyperthreading in BLAS
doFuture::registerDoFuture()
future::plan(future::multisession, workers = n_threads)
print(paste0("Running sampling test 1 with ", future::nbrOfWorkers(), " workers."))


# run simulation ----------------------------------------------------------
settings = expand.grid(method = methods, dimension = dimensions)
n_settings = nrow(settings)

print( paste0("*** Total methods: ", n_settings, " ***"))

progressr::with_progress({
  p = progressr::progressor(along = 1:(n_settings * n_dims))
  
  foreach(i = 1:n_settings, .inorder = FALSE, 
          .export = ls(globalenv())) %dorng% {
            
    p(sprintf("i=%g", i))
    
    method = settings[[i, "method"]]
    d = settings[[i, "dimension"]]
    
    registerDoRNG(seed)
    
    # problem-dimension specific settings
    problem_params = list(
      n = 1000,
      mu = rep(0, d),
      Sigma = diag(d),
      lb = rep(-2, d),
      ub = rep(2, d),
      initial = rep(0, d)
    )
    
    # method specific settings
    param_string = paste(
      names(method$parameters),
      method$parameters, 
      sep = "=", collapse = ", "
    )
    
    print( paste0("--- ", method$method, ", Dimension: ", d, "---") )
    print( paste("Parameters: ", param_string) )
    
    # all input arguments
    params = c(problem_params, method$parameters)
    
    
    tictoc::tic()
    samples = do.call(method$method, params)
    elapsed = tictoc::toc(quiet = TRUE)
        
    runtime = data.frame(
      method = method$method,
      runtime = elapsed$toc - elapsed$tic,
      n_dim = d
    )
    
    # save samples
    saveRDS(samples, paste0(
      sample_path, "/", method$method, "_d=", d
    ))
    
    # save timings
    saveRDS(runtime, paste0(
      timing_path, "/", method$method, "_d=", d
    ))
  } 
})



