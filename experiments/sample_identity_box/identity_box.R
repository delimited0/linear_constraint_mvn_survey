'identity_box

Usage:
  identity_box.R (--method_conf=<method_conf>) (--dim_conf=<dim_conf>) (--half_width=<half_width>) (--sample_path=<sample_path>) [--seed=seed] [--n_threads=n_threads] 
  identity_box.R (-h|--help)

Options:
  -h --help  Usage.
  --method_conf=<method_conf>  Method configuration.
  --dim_conf=<dim_conf>  Dimension configuration.
  --half_width=<half_width>  Truncation box half width.
  --sample_path=<sample_path>  Configuration specific sampler output.
  --seed=seed  Seed.
  --n_threads=n_threads  Number of cores.
' -> doc

opts = docopt::docopt(doc, version = 'identity_box 1.0')

method_conf = opts$method_conf
dim_conf = opts$dim_conf
half_width = as.numeric(opts$half_width)
sample_path = opts$sample_path
seed = as.numeric(opts$seed)
n_threads = as.numeric(opts$n_threads)

# default arguments -------------------------------------------------------
if (is.null(seed)) seed = 2022
if (is.null(n_threads)) n_threads = 1

# hard coded arguments for debugging --------------------------------------
# method_conf = "experiments/sample_identity_box/method_conf.json"
# dim_conf = "experiments/sample_identity_box/test_dim_conf.json"
# sample_path = "experiments/sample_identity_box/samples"
# half_width = 3

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
future::plan(future::multicore, workers = n_threads)


# start up summary --------------------------------------------------------
print(paste0("Running identity box sampling with ",
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
  
  foreach(i = 1:n_settings, .inorder = FALSE, .options.RNG = seed,
          .export = ls(globalenv())) %dorng% 
    {
      
      method = settings[[i, "method"]]
      d = settings[[i, "dimension"]]
      
      p(message = sprintf("%s, dimension %d", method, d))
      
      # problem-dimension specific settings
      problem_params = list(
        n = 5000,
        mu = rep(0, d),
        Sigma = diag(d),
        lb = rep(-half_width, d),
        ub = rep(half_width, d),
        initial = rep(0, d)
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
      attr(samples, "half_width") = half_width
      
      # handle output directories
      method_sample_path = 
        paste0(sample_path, "/width=", half_width, "/", method$method, "/")
      if (!dir.exists(method_sample_path)) 
        dir.create(method_sample_path, recursive=TRUE)
      
      # save samples
      saveRDS(samples, paste0(method_sample_path, "d=", d))
    } 
}, enable = TRUE)


