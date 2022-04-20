'proton_radius

Usage:
  proton_radius.R (--method_conf=<method_conf>) (--sample_path=<sample_path>) (--reps=<reps>) [--seed=seed] [--n_threads=n_threads] 
  proton_radius.R (-h|--help)

Options:
  -h --help  Usage.
  --method_conf=<method_conf>  Method configuration.
  --sample_path=<sample_path>  Configuration specific sampler output.
  --reps=<reps>  Number of repetitions.
  --seed=seed  Seed.
  --n_threads=n_threads  Number of cores.
' -> doc

opts = docopt::docopt(doc, version = 'proton_radius 1.0')

method_conf = opts$method_conf
sample_path = opts$sample_path
# timing_path = opts$timing_path
seed = as.numeric(opts$seed)
n_threads = as.numeric(opts$n_threads)
n_reps = as.numeric(opts$reps)

# default arguments -------------------------------------------------------
if (is.null(seed)) seed = 2022
if (is.null(n_threads)) n_threads = 1

# hard coded arguments for debugging --------------------------------------
# method_conf = "experiments/proton_radius/method_conf.json"

# libraries ---------------------------------------------------------------
library(here)
library(doFuture)
library(progressr)
library(doRNG)
library(foreach)
library(tictoc)

source(here("sampling_wrapper.R"))
source(here("experiments", "proton_radius", "paper_sampler", "cGP.R"))
source(here("experiments", "proton_radius", "paper_sampler", "maternCov.R"))


# read settings ---------------------------------------------------------
methods = jsonlite::read_json(method_conf, simplifyVector = FALSE)

# setup output directories -------------------------------------------
if (!dir.exists(sample_path)) dir.create(sample_path)

# parallel set up ---------------------------------------------------------
RhpcBLASctl::blas_set_num_threads(1)  # no hyperthreading in BLAS
doFuture::registerDoFuture()
future::plan(future::multisession, workers = n_threads)

# simulation setting ------------------------------------------------------

# basis functions
h = function(x, j)
  ifelse(abs(x-u[j]) <= dN, (1 - abs(x - u[j]) / dN), 0)

phi = function(x,j) {
  if(x < u[j]-dN) { 
    0
  }
  else if (u[j]-dN<=x & x<u[j]) {
    (x+dN-u[j])^2/2/dN
  } 
  else if (u[j]<=x & x<u[j]+dN) {
    dN - (x-dN-u[j])^2/2/dN
  } else{
    dN
  }
} 

phi_1 = function(x) 
  ifelse(x <= u[1]+dN, dN/2 - (dN+u[1]-x)^2/(2*dN), dN/2)
phi_N = function(x) 
  ifelse(x < u[N]-dN, 0, (x+dN-u[N])^2/(2*dN))

psi = function(x,j) {
  if(x < u[j]-dN) {
    0
  }
  else if (u[j]-dN<=x & x<u[j]) {
    (x+dN-u[j])^3/(6*dN)
  }
  else if (u[j]<=x & x<u[j]+dN) {
    dN*(x-u[j]) - (x-dN-u[j])^3/(6*dN)
  }
  else {
    dN^2 + dN*(x-dN-u[j]) 
  }
}

psi_1 = function(x) 
  ifelse(x < dN & x!=0, dN*x/2 - (x-dN)^3/(6*dN), dN^2/2 + dN*(x-dN)/2)
psi_N = function(x) 
  ifelse(x < u[N]-dN, 0, (x+dN-u[N])^3/(6*dN))

# read real experimental covariates
Q2 = read.table(here("experiments", "proton_radius", "xvals.txt"))$V1
Q2max = max(Q2)

# dimensionless scaled variable
x = Q2 / Q2max  
n = length(x)

# dipole function
rp = 0.84  # ground truth radius
p1 = 12.0 / rp^2
y_true = (1 + (x / p1))^(-2)

# Gaussian error standard deviation
sig = 0.005

# basis matrix
C = max(x)
# number of knots 
N = floor(n/4)+1 
# define equal-spaced knots 
u = seq(0,C, length.out = N)
dN = C/(N-1)
Phi_x = matrix(nrow = n, ncol = N)
for(i in 1:n){
  Phi_x[i,1] = psi_1(x[i])
  Phi_x[i,N] = psi_N(x[i])
  for(j in 2:(N-1)){Phi_x[i,j] = psi(x[i],j)}
}
Phi = cbind(as.matrix(rep(1,n),nrow=n), x, Phi_x)

# transformation coefficient matrices
max_phi = c(dN/2, rep(dN, N-2), dN/2)
# for cGP, c0GP and uGP
trans_mat = matrix(nrow = N+1, ncol = N+1) 
trans_mat[1,] = -c(1, max_phi)
trans_mat[-1,] = cbind(as.matrix(rep(0,N),ncol=1),diag(N)) 

# smoothness parameter 
nu = 2.5
# scale-length parameter
l = 20
# number of mcmc iterations
Niter = 500

# start up summary --------------------------------------------------------
n_methods = length(methods)

print(paste0("Running proton radius estimation with ",
             future::nbrOfWorkers(), " workers."))
print(paste0("Comparing ", n_methods, " methods"))

cat( paste0( sapply(methods, function(x) x$method) , collapse = "\n") )
cat("\n")

# run simulation ----------------------------------------------------------

for (i in 1:n_methods) {
  
  method = methods[[i]]$method
  params = methods[[i]]$parameters
  
  method_sample_path = paste0(sample_path, "/", method, "/")
  if (!dir.exists(method_sample_path)) dir.create(method_sample_path)
  
  print( paste0("--- ", method, "---") )
    
  registerDoRNG(seed)
  progressr::with_progress({
    p = progressr::progressor(along = 1:(n_methods))    
    
    foreach(i = 1:n_reps, .inorder = FALSE, 
            .export = ls(globalenv()), .errorhandling = "remove") %dorng% {
              
      p(sprintf("i=%g", i))
              
      # simulate data
      y_obs = y_true + rnorm(n, 0, sig)
      
      tic()
      posterior = cGP(x, y_obs,  nu = nu,  l = l*C,
                   niter = Niter, u = u, Phi = Phi, trans_mat = trans_mat,
                   sampler = method, sampler_params = params)
      elapsed = toc(quiet=TRUE)
      
      samples = t(posterior$weights)
      
      attr(samples, "method") = method
      attr(samples, "runtime") = elapsed$toc - elapsed$tic
      
      saveRDS(samples, paste0(method_sample_path, "rep=", i))
    }
  })
}
  