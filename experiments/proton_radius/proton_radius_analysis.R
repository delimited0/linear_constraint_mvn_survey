library(here)
library(data.table)
library(coda)
library(mcmcse)

Q2 = read.table(here("experiments", "proton_radius", "xvals.txt"))$V1
Q2max = max(Q2)
x = Q2 / Q2max
n = length(x)

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

# test samples ------------------------------------------------------------

sample_path = here("experiments", "proton_radius", "test_samples")

all_stats = rbindlist(lapply(
  dir(sample_path, full.names = TRUE), 
  function(method_path) {
    
    rbindlist(lapply(
      dir(method_path, full.names = TRUE), 
      function(method_rep) {
        
        samples = readRDS(method_rep)
        d = ncol(samples)
        
        # how does Q_max enter here?
        rp_samples = sqrt(-6 * samples[-c(1:100), 2]) 
        
        point_est_summary = mcse(rp_samples)
        
        data.table(
          posterior_mean = point_est_summary$est,
          posterior_median = median(rp_samples),
          mc_se = point_est_summary$se,
          runtime = attr(samples, "runtime"),
          method = attr(samples, "method"),
          ess = effectiveSize(rp_samples)
        )
      }
    ))     
  }
))

# visualization -----------------------------------------------------------
library(ggplot2)

ggplot(all_stats, aes(x = method, y = posterior_median)) +
  geom_point()

ggplot(all_stats, aes(x = method, y = mc_se)) +
  geom_point()

ggplot(all_stats, aes(x = method, y = ess))

