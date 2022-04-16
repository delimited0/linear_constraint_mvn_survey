rm(list = ls())

library(MASS)
library(FastGP)
library(mcmcplots)
library(fields)
library(HDInterval)
library(mcmcse)
library(ggplot2)

source('RobustGP.R')
source('ESS_pro.R')
source('ESS_joint.R')
source('function_pro.R')

#####################
# Inputs parameters #
#####################
# Input data 
Y = read.table("YF0.txt")
X = read.table("xvals.txt")

# sample size
n = nrow(X)

# standard deviation of normal errors for Y
sig.pro = 0.001
# smoothness parameter for Matern kernel
nu = 0.5
# number of basis function
N = ceiling(n/8)-1
# length-scale parameter
l = 10 
# hyperparameter of the normalizing factor \xi_0
sig0 = 0.0005
# specifying the probability within the credible interval
prob = 0.95
# seeds
sseed = 67
# number of MCMC iterations
mcmc = 10000 
# number of burn-in samples
brn = 4000
# thining 
thin = 10

#########################
# Running the algorithm #
#########################

#c_1GP
Results = RobustGP(Y,X,sig=sig.pro,nu=nu,N=N,l=l,sig0=sig0,prob=prob,sseed=sseed,
                   mcmc=mcmc,brn=brn,thin=thin,return.plot = TRUE,return.traplot = FALSE)

#cGP (set xi0.fix=1)
#Results = RobustGP(Y,X,sig=sig.pro,nu=nu,N=N,l=l,sig0=sig0,prob=prob,sseed=sseed,
#                   mcmc=mcmc,brn=brn,thin=thin,xi0.fix=1,return.plot = TRUE,return.traplot = FALSE)



###########
# Outputs #
###########

#Results: a list object including the following:
  
# (1) Estimated Radius: Posterior median:
r_est = Results[[1]] 


#(2) $95\%$ confidence interval of the proton radius:
r_CI = Results[[2]]


# (3) MCMC samples of the proton radius:
r_sam = Results[[3]]


# (4) Data of $G_E$ and $Q^2$, saved for replications 
data = Results[[4]]


# (5) Posterior samples on $\xi,\xi_0,\xi_1,\tau^2,\sigma^2$ and posterior mean as well as $95\%$ CI of $G_E$ estimates:
post.sam = Results[[5]]


# (6) Watanabe's AIC (a model choice criterion) for the current choice of hyperparameters:
model.check = Results[[6]]
WAIC = model.check[[1]][1]


# (7) The shortest set (interval in this case) with "prob"% posterior probability:
HD_rp = Results[[7]]


# (8) Estimate effective sample size for the MCMC path: 
r_Ess = Results[[8]]


# (9) Monte Carlo standard error:
r_mcse = Results[[9]]
# The estimate of expectation of proton radius:
r_estimate = r_mcse[1]
# Standard error of MCMC samples 
r_error = r_mcse[2]
# Standard error of data, MCMC iteration stops when the Monte Carlo standard error is small compared to the variability in the target distribution
sd(as.numeric(as.matrix(Y)))



#END
