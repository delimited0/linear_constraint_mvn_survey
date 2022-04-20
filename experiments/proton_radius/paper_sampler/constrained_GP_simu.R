library(MASS)
library(Matrix)
library(mvtnorm)
library(invgamma)
library(TruncatedNormal)
library(ggplot2)
library(grid)

source(here::here("experiments", "proton_radius", "paper_sampler", "maternCov.R")) # matern kernel with gerneral smoothness \nu and length-scale parameter l
source(here::here("experiments", "proton_radius", "paper_sampler", "cGP.R"))
source(here::here("experiments", "proton_radius", "paper_sampler", "c0GP.R"))
source(here::here("experiments", "proton_radius", "paper_sampler", "c1GP.R"))
source(here::here("experiments", "proton_radius", "paper_sampler", "uGP.R"))

source(here::here("sampling_wrapper.R"))

## basis functions ##
h = function(x,j) ifelse(abs(x-u[j]) <= dN, (1-abs(x-u[j])/dN), 0)

phi = function(x,j) if(x < u[j]-dN){ 0
}else if(u[j]-dN<=x & x<u[j]){ (x+dN-u[j])^2/2/dN
}else if(u[j]<=x & x<u[j]+dN){ dN - (x-dN-u[j])^2/2/dN
}else{dN}

phi_1 = function(x) ifelse(x <= u[1]+dN, dN/2 - (dN+u[1]-x)^2/(2*dN), dN/2)
phi_N = function(x) ifelse(x < u[N]-dN, 0, (x+dN-u[N])^2/(2*dN))

psi = function(x,j) if(x < u[j]-dN){ 0
}else if(u[j]-dN<=x & x<u[j]){ (x+dN-u[j])^3/(6*dN)
}else if(u[j]<=x & x<u[j]+dN){dN*(x-u[j]) - (x-dN-u[j])^3/(6*dN)
}else{ dN^2 + dN*(x-dN-u[j]) }

psi_1 = function(x) ifelse(x < dN & x!=0, dN*x/2 - (x-dN)^3/(6*dN), dN^2/2 + dN*(x-dN)/2)
psi_N = function(x) ifelse(x < u[N]-dN, 0, (x+dN-u[N])^3/(6*dN))



## simulate data -- (x, y) ##
a1 = -0.1176 
D = function(x) 1/(1-a1*x/2)^2 
n = 250
x = runif(n, min = 0, max = 10)
sigma = 0.002
err = sigma*rnorm(n,0,1)
y = D(x) + err

## basis matrix ##
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

## transformation coefficient matrices ##
max_phi = c(dN/2, rep(dN, N-2), dN/2)
# for cGP, c0GP and uGP
trans_mat = matrix(nrow = N+1, ncol = N+1) 
trans_mat[1,] = -c(1, max_phi)
trans_mat[-1,] = cbind(as.matrix(rep(0,N),ncol=1),diag(N)) 
# for c1GP
trans_mat1 = as.matrix(bdiag(1, trans_mat))

# smoothness parameter 
nu = 2.5
# scale-length parameter
l = 15
# number of mcmc iterations
Niter = 500


# cGP #
r = cGP(x, y,  nu = nu,  l = l*C,
        niter = Niter, u = u, Phi = Phi, trans_mat = trans_mat,
        sampler = epess, sampler_params = list(n_perthresh = 1, n_slice = 1))

r = cGP(x, y,  nu = nu,  l = l*C,
        niter = Niter, u = u, Phi = Phi, trans_mat = trans_mat,
        sampler = gibbs_lg2015, sampler_params = list(tuvn_sampler = "lg2015"))

r = cGP(x, y,  nu = nu,  l = l*C,
        niter = Niter, u = u, Phi = Phi, trans_mat = trans_mat,
        sampler = rsm)

r = cGP(x, y,  nu = nu,  l = l*C,
        niter = Niter, u = u, Phi = Phi, trans_mat = trans_mat,
        sampler = gibbs_cov)

r = cGP(x, y,  nu = nu,  l = l*C,
        niter = Niter, u = u, Phi = Phi, trans_mat = trans_mat,
        sampler = gibbs_prec)

r = cGP(x, y,  nu = nu,  l = l*C,
        niter = Niter, u = u, Phi = Phi, trans_mat = trans_mat,
        sampler = rhmc, sampler_params = list(traj_length = 2))

# c0GP #
r0 = c0GP(x, y,  nu = nu , l = l*C,
          niter = Niter, u = u, Phi = Phi, trans_mat = trans_mat)

# c1GP #
r1 = c1GP(x, y,  nu = nu , l = l*C,
          niter = Niter, u = u, Phi = Phi, trans_mat = trans_mat1)

r1 = c1GP(x, y,  nu = nu , l = l*C,
          niter = Niter, u = u, Phi = Phi, trans_mat = trans_mat1,
          sampler = rsm)

r1 = c1GP(x, y,  nu = nu , l = l*C,
          niter = Niter, u = u, Phi = Phi, trans_mat = trans_mat1,
          sampler = gibbs_cov)

r1 = c1GP(x, y,  nu = nu , l = l*C,
          niter = Niter, u = u, Phi = Phi, trans_mat = trans_mat1,
          sampler = gibbs_prec)

r1 = c1GP(x, y,  nu = nu , l = l*C,
          niter = Niter, u = u, Phi = Phi, trans_mat = trans_mat1,
          sampler = rhmc, sampler_params = list(traj_length = 2))

r1 = c1GP(x, y,  nu = nu , l = l*C,
          niter = Niter, u = u, Phi = Phi, trans_mat = trans_mat1,
          sampler = ghk)



# uGP #
ru = constrGP_n(x, y,  nu = nu , l = l*C,
                niter = Niter, u = u, Phi = Phi, trans_mat = trans_mat)


## posterior median estimates and 95% posterior credible intervals ## 

f = r[[1]]
f_pred = apply(Phi%*%f[,-c(1:100)],1, median)
f_CI = apply(Phi%*%f[,-c(1:100)],1, quantile, probs = c(0.025,0.975))

# cGP #
f = r0[[1]]
f_pred = apply(Phi%*%f[,-c(1:100)],1, median)
f_CI = apply(Phi%*%f[,-c(1:100)],1, quantile, probs = c(0.025,0.975))
# c0GP #
fn = r0_n[[1]]
fn_pred = apply(Phi%*%fn[,-c(1:100)],1, median)
fn_CI = apply(Phi%*%fn[,-c(1:100)],1, quantile, probs = c(0.025,0.975))
# c1GP  #
fn1 = r_n1[[1]]
fn1_pred = apply(Phi%*%fn1[,-c(1:100)],1, median)
fn1_CI = apply(Phi%*%fn1[,-c(1:100)],1, quantile, probs = c(0.025,0.975))
# uGP # 
fu = ru_n[[1]]
fu_pred = apply(Phi%*%fu[,-c(1:100)],1, median)
fu_CI = apply(Phi%*%fu[,-c(1:100)],1, quantile, probs = c(0.025,0.975))



## plot on test grid ##
xtest = seq(0, max(x), length.out = n)
ytest = D(xtest)

ntest = length(xtest)
Phi_test = matrix(nrow = ntest, ncol = N)
for(i in 1:ntest){
  Phi_test[i,1] = psi_1(xtest[i])
  Phi_test[i,N] = psi_N(xtest[i])
  for(j in 2:(N-1)){Phi_test[i,j] = psi(xtest[i],j)}
}

Phi_t = cbind(as.matrix(rep(1,ntest),nrow=ntest), xtest, Phi_test)
f_test = apply(Phi_t%*%f[,-c(1:100)],1, median)
CI_test = apply(Phi_t%*%f[,-c(1:100)],1, quantile, probs = c(0.025,0.975))

fn_test = apply(Phi_t%*%fn[,-c(1:100)],1, median)
CIn_test = apply(Phi_t%*%fn[,-c(1:100)],1, quantile, probs = c(0.025,0.975))

fn1_test = apply(Phi_t%*%fn1[,-c(1:100)],1, median)
CIn1_test = apply(Phi_t%*%fn1[,-c(1:100)],1, quantile, probs = c(0.025,0.975))

fu_test = apply(Phi_t%*%fu[,-c(1:100)],1, median)
CIu_test = apply(Phi_t%*%fu[,-c(1:100)],1, quantile, probs = c(0.025,0.975))

df = data.frame(xtest,f_test)
G <- ggplot(df, aes(x=xtest, y=f_test), color=variable)+ 
  geom_line(aes(x=xtest, y=f_test), colour="blue") 
+
  # geom_line(aes(x=xtest, y=fn_test), colour="green")+
  # geom_line(aes(x=xtest, y=fu_test), colour="purple")+
  # geom_line(aes(x=xtest, y=ytest), colour="red", lty = 2)+
  # geom_line(aes(x=xtest, y=CIu_test[1,]), colour="purple", lty = 2)+
  # geom_line(aes(x=xtest, y=CIu_test[2,]), colour="purple", lty = 2)+ 
  # geom_line(aes(x=xtest, y=CIn1_test[1,]), colour="blue", lty = 5, cex = 0.3)+
  # geom_line(aes(x=xtest, y=CIn1_test[2,]), colour="blue", lty = 5, cex = 0.3)+ 
  # geom_ribbon(aes(ymin=CI_test[1,], ymax=CI_test[2,]), alpha=0.4)+
  # geom_ribbon(aes(ymin=CIn_test[1,], ymax=CIn_test[2,]), alpha=0.2)

G <- G + theme_bw()
Glabs <- G+labs(x = "x", y = "y")
Glabs + theme( 
  axis.title.x = element_text(size=8, face="bold"),
  axis.title.y = element_text(size=8, face="bold"),
  legend.text = element_text(colour="blue", size=10, face="bold"))

Gz <- G + coord_cartesian(xlim = c(0,0.1), ylim = c(0.985,1.0), expand = TRUE)
Gzlabs<- Gz+labs(x = "x", y = "y")
Gzlabs2<- Gzlabs + theme( 
  axis.title.x = element_text(size=8, face="bold"),
  axis.title.y = element_text(size=8, face="bold"),
  legend.text = element_text(colour="blue", size=10, face="bold"))

vp <- viewport(width = 0.45, height = 0.45, x = 0.75, y = 0.75)
print(Gzlabs2, vp=vp)

rp_samples = sqrt(-6* r$weights[2,-c(1:100)])
point_est_summary = mcse(rp_samples)
