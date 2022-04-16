#########################################################
# GP model incorporating constraints (1).
#
# Input: x:         Covariate, a n*1 vector. 
#        y:         Response, a n*1 vector.  
#        nu:        Smoothness parameter in matern kernel. 
#        l:         Length-scale parameter in matern kernel.
#        niter:     Iteration numvers of MC. 
#        u:         Knots of basis functions.
#        Phi:       Basis function evaluation, a n*(N+2) matrix.
#        trans_mat: tansformation matrix used in the constrained GP methods.
#
# Output: f:        mcmc samples of weights parameter.  
#         tau:      mcmc samples of signal-to-noise level in matern kernel. 
#         sigma:    mcmc samples of noise level.
#
##########################################################


c0GP = function(x, y, nu, l, niter, u, Phi, trans_mat){
  
  n=length(x)
  
  # covariance matrix #  
  mcov = maternCov(u, l, nu = nu)
  D = eigen(mcov)$values
  U = eigen(mcov)$vectors
  for(i in 1:length(D)){ if(D[i] < 0) D[i] = -D[i]}
  mcov = U%*%diag(D)%*%t(U)
  
  # initial values #
  tau = c()
  tau[1] = 0.01 
  sigma = c()
  sigma[1] = 0.01
  mcov_inv = chol2inv(chol(mcov))
  f = matrix(nrow = N+2, ncol = niter)
  f[,1] = t(chol(mcov))%*%rnorm(N+2, 0, 1)*sqrt(tau[1])
  re = list()
  
  # mcmc runs #
  for(i in 2:niter){
    
    if(i%%10 == 0) print(i)
    
    # update f #
    Sig = tau[i-1]*t(Phi) %*% Phi + sigma[i-1]*mcov_inv
    Siginv = tau[i-1]*sigma[i-1]*chol2inv(chol(Sig))
    postmean = Siginv%*%t(Phi)%*%y/sigma[i-1]
    mu1 = postmean[1]
    mu2 = postmean[-1]
    Siginv_11 = Siginv[1,1]
    Siginv_21 = as.matrix(Siginv[-1,1])
    Siginv_12 = t(Siginv_21) 
    Siginv_22 = Siginv[-1,-1]
    
    mu_cond = mu2 + Siginv_21/Siginv_11*(1-mu1)
    Sig_cond = (Siginv_22 - Siginv_21%*%Siginv_12/Siginv_11)
    f_cond = mu_cond + t(chol(Sig_cond))%*%rnorm(N+1,0, 1)
    f[,i] = as.matrix(c(1,f_cond), ncol = 1)
    
    # update sigma #
    sigma[i] = rinvgamma(1, shape = n/2, rate = t(y-Phi%*%f[,i])%*%(y-Phi%*%f[,i])/2)
    
    # update tau #
    tau[i] = rinvgamma(1, shape = (N+2)/2, rate = t(f[,i])%*%mcov_inv%*%f[,i]/2)
  }
  
  re = list(f, tau, sigma)
  
  return(re)
  
}

