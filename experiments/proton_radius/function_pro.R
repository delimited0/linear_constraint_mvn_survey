### All required functions for using ESS
### Functions related to Wood and Chan algorithm of drawing samples
### MH for sampling from \nu and \ell
### Covariance matrix and design matrix (using basis function) are also defined
### And all related and dependant functions are here

### Required libraries:
library(fields)
library(FastGP)

### Define design matrix ###

### The basis functions ###

# For monotone function estimation:
psi_j=function(x,my_knot,delta_N)
{
  N=length(my_knot)
  k=rep(0,N)
  i=max(which(my_knot<=x))
  if(i==1)
  {
    k[1]=x-0.5*(x^2)/delta_N
    k[2]=x-my_knot[2]*x/delta_N+0.5*x^2/delta_N
  }
  if(i==2)
  {
    k[1]=delta_N/2
    k[2]=delta_N/2+(x-my_knot[2])*(1+my_knot[2]/delta_N)-0.5*(x^2-my_knot[2]^2)/delta_N
    k[3]=(x-my_knot[2])*(1-my_knot[3]/delta_N)+0.5*(x^2-my_knot[2]^2)/delta_N
  }
  if(i==N)
  {
    k[1]=delta_N/2
    k[2:(N-1)]=delta_N
    k[N]=delta_N/2
  }
  if(i!=1 && i!=2 && i!=N)
  {
    k[1]=delta_N/2
    k[2:(i-1)]=delta_N
    k[i]=delta_N/2+(x-my_knot[i])*(1+my_knot[i]/delta_N)-0.5*(x^2-my_knot[i]^2)/delta_N
    k[i+1]=(x-my_knot[i])*(1-my_knot[i+1]/delta_N)+0.5*(x^2-my_knot[i]^2)/delta_N
  }
  return(k)
}

# For convex function estimation:
phi_j=function(x,my_knot,delta_N)
{
  N=length(my_knot)
  k=rep(0,N)
  if(x>=my_knot[1] && x<my_knot[2])
  {
    k[1]=(x^2/2)-0.5*(x^3/3)/delta_N
    k[2]=0.5*(x^3/3)/delta_N
  }
  if(x>=my_knot[2] && x<my_knot[3])
  {
    k[1]=(my_knot[2]^2/2)-0.5*(my_knot[2]^3/3)/delta_N+0.5*delta_N*(x-my_knot[2])
    k[2]=my_knot[2]^3/(3*delta_N)+0.5*delta_N*(x-my_knot[2])+(x^2-my_knot[2]^2)-1.5*my_knot[2]*(x-my_knot[2])-0.5*(x^3/3)/delta_N
    k[3]=(1-my_knot[3]/delta_N)*(x^2/2-x*my_knot[2]+my_knot[2]^2/2)+0.5*(x^3/3-x*my_knot[2]^2+2*my_knot[2]^3/3)/delta_N
  }
  if(x>=my_knot[3] && x<my_knot[N])
  {
    k[1]=(my_knot[2]^2/2)-0.5*(my_knot[2]^3/3)/delta_N+0.5*delta_N*(x-my_knot[2])
    k[2]=my_knot[2]^3/(3*delta_N)+0.5*delta_N^2+my_knot[3]^2-2.5*my_knot[2]^2-0.5*(my_knot[3]^3/3)/delta_N+delta_N*(x-my_knot[3])
    k[3]=(1-my_knot[3]/delta_N)*(my_knot[3]^2/2-my_knot[3]*my_knot[2]+my_knot[2]^2/2)+0.5*(my_knot[3]^3/3-my_knot[3]*my_knot[2]^2+2*my_knot[2]^3/3)/delta_N+delta_N*(x-my_knot[3])
    
    if(N>4){
      for(j in 4:(N-1)){
        if(x<my_knot[j-1])
          k[j]=0
        if(x>=my_knot[j-1] && x<my_knot[j])
          k[j]=(1-my_knot[j]/delta_N)*(0.5*(x^2-my_knot[j-1]^2)-my_knot[j-1]*(x-my_knot[j-1]))+0.5*((x^3/3-my_knot[j-1]^3/3)-my_knot[j-1]^2*(x-my_knot[j-1]))/delta_N
        if(x>=my_knot[j] && x<my_knot[j+1])
          k[j]=(1-my_knot[j]/delta_N)*(0.5*(my_knot[j]^2-my_knot[j-1]^2)-my_knot[j-1]*(my_knot[j]-my_knot[j-1]))+0.5*((my_knot[j]^3/3-my_knot[j-1]^3/3)-my_knot[j-1]^2*(my_knot[j]-my_knot[j-1]))/delta_N+0.5*delta_N*(x-my_knot[j])+
            (1+my_knot[j]/delta_N)*(0.5*(x^2-my_knot[j]^2)-my_knot[j]*(x-my_knot[j]))-0.5*((x^3/3-my_knot[j]^3/3)-my_knot[j]^2*(x-my_knot[j]))/delta_N
        if(x>=my_knot[j+1])
          k[j]=(1-my_knot[j]/delta_N)*(0.5*(my_knot[j]^2-my_knot[j-1]^2)-my_knot[j-1]*(my_knot[j]-my_knot[j-1]))+0.5*((my_knot[j]^3/3-my_knot[j-1]^3/3)-my_knot[j-1]^2*(my_knot[j]-my_knot[j-1]))/delta_N+0.5*delta_N*(my_knot[j+1]-my_knot[j])+
            (1+my_knot[j]/delta_N)*(0.5*(my_knot[j+1]^2-my_knot[j]^2)-my_knot[j]*(my_knot[j+1]-my_knot[j]))-0.5*((my_knot[j+1]^3/3-my_knot[j]^3/3)-my_knot[j]^2*(my_knot[j+1]-my_knot[j]))/delta_N+delta_N*(x-my_knot[j+1])
      }
    }
    
    if(x>=my_knot[N-1] && x<my_knot[N])
      k[N]=(1-my_knot[N]/delta_N)*(0.5*(x^2-my_knot[N-1]^2)-my_knot[N-1]*(x-my_knot[N-1]))+0.5*((x^3/3-my_knot[N-1]^3/3)-my_knot[N-1]^2*(x-my_knot[N-1]))/delta_N
  }
  if(x>=my_knot[N])
  {
    k[1]=(my_knot[2]^2/2)-0.5*(my_knot[2]^3/3)/delta_N+0.5*delta_N*(x-my_knot[2])
    k[2]=my_knot[2]^3/(3*delta_N)+0.5*delta_N^2+my_knot[3]^2-2.5*my_knot[2]^2-0.5*(my_knot[3]^3/3)/delta_N+delta_N*(x-my_knot[3])
    k[3]=(1-my_knot[3]/delta_N)*(my_knot[3]^2/2-my_knot[3]*my_knot[2]+my_knot[2]^2/2)+0.5*(my_knot[3]^3/3-my_knot[3]*my_knot[2]^2+2*my_knot[2]^3/3)/delta_N+delta_N*(x-my_knot[3])
    for(j in 4:(N-1)){
      k[j]=(1-my_knot[j]/delta_N)*(0.5*(my_knot[j]^2-my_knot[j-1]^2)-my_knot[j-1]*(my_knot[j]-my_knot[j-1]))+0.5*((my_knot[j]^3/3-my_knot[j-1]^3/3)-my_knot[j-1]^2*(my_knot[j]-my_knot[j-1]))/delta_N+0.5*delta_N*(my_knot[j+1]-my_knot[j])+
        (1+my_knot[j]/delta_N)*(0.5*(my_knot[j+1]^2-my_knot[j]^2)-my_knot[j]*(my_knot[j+1]-my_knot[j]))-0.5*((my_knot[j+1]^3/3-my_knot[j]^3/3)-my_knot[j]^2*(my_knot[j+1]-my_knot[j]))/delta_N+delta_N*(x-my_knot[j+1])
    }
    k[N]=(1-my_knot[N]/delta_N)*(0.5*(my_knot[N]^2-my_knot[N-1]^2)-my_knot[N-1]*(my_knot[N]-my_knot[N-1]))+0.5*((my_knot[N]^3/3-my_knot[N-1]^3/3)-my_knot[N-1]^2*(my_knot[N]-my_knot[N-1]))/delta_N+0.5*delta_N*(x-my_knot[N])
  }
  return(k)
}

### Function to form design matrix:
des.mat1=function(x,my_knot,delta_N){
  # Function to form basis matrix for monotone constraint
  n=length(x)
  N=length(my_knot)-1
  # design matrix \Psi(n X N+1)
  X=matrix(0,n,N+1)
  for(l in 1:n){
    X[l,1:(N+1)]=psi_j(x[l],my_knot,delta_N)
  }
  return(X)
}

des.mat2=function(x,my_knot,delta_N){
  # Function to form basis matrix for convex constraint
  n=length(x)
  N=length(my_knot)-1
  # design matrix \Phi(n X N+1)
  X=matrix(0,n,N+1)
  for(l in 1:n){
    X[l,1:(N+1)]=phi_j(x[l],my_knot,delta_N)
  }
  return(X)
}

### Function to compute \xi_1 * x :
xix=function(a,b){
  return(a*b)
}

#Given a \nu (smoothness parameter of matern kernel) finding a value of 
# l (length-scale parameter) such that the correlation between the 
# maximum seperation is some small value, say 0.05

# Matern kernel with smoothness nu and length-scale l:
MK = function(x, y ,l, nu){
  ifelse(abs(x-y)>0, (sqrt(2*nu)*abs(x-y)/l)^nu/(2^(nu-1)*gamma(nu))*besselK(x=abs(x-y)*sqrt(2*nu)/l, nu=nu), 1.0)
}

# Integration of Matern kernel with nu and l (assume x \in [-a,1]):
MK_int_x = function(x, y, l, nu, nt, a){
  delta_t = seq(-a, x, length.out = nt)
  val = sum(MK(delta_t, y, l, nu))*x/nt
  return(val)     
}

# Double integrartion of Matern kernel with nu and l (assume x \in [-a,1]):
MK_int_xy = function(x, y, l, nu, nt, a){
  delta_t1 = seq(-a, x, length.out = nt)
  delta_t2 = seq(-a, y, length.out = nt)
  val=c()
  for(i in 1:length(delta_t2)){
   val[i] = sum(MK(delta_t1, delta_t2[i], l, nu))*(x+a)/nt*(y+a)/nt}
   Val = sum(val)
  return(Val)     
}

# function for uniroot:
fl=function(l,para){ 
  #para[1]=x, para[2]=y and para[3]=nu of MK : Matern kernel function;
  #para[4]=pre-specified value of the correlation
  a=MK(para[1],para[2],l,para[3])
  return(a-para[4])
}

# function for estimating l:
l_est=function(nu,range,val){
  # nu : smoothness; range : c(min, max) of the range of variable
  # val : pre-specified value of the correlation between the maximum seperation
  para=c(range[1],range[2],nu,val)
  rl=uniroot(f=fl,interval = c(0.000001,100000),para)
  return(rl$root)
}

# Covariance matrix
covmat=function(knot,nu,l){
  return(MK(rdist(knot),0,l,nu))
}

# Covariance mateix for the joint distribution
# nt : number of partitions for numerical calulation
#a : boundary for numerical calulation

covmat_joint = function(knot,l,nu,nt,a){
  nknot = length(knot)
  covmat1 = covmat(knot,nu,l)
  Covmat = matrix(0,nrow=(nknot+1), ncol=(nknot+1))
  Covmat[1,1] = MK_int_xy(0,0,l,nu,nt,a)
  for(i in 1:nknot){Covmat[1,(i+1)] = MK_int_x(0,knot[i],l,nu,nt,a)}
  Covmat[-1,1] = t(Covmat[1,-1])
  Covmat[-1,-1] = covmat1
  return(Covmat)
}

# Order of the circulant matrix:
# minimum value of g and m so that G can be embedded into C
min_g=function(knot){
  N=length(knot)
  g=ceiling(log(2*N,2))   #m=2^g and m>=2(n-1) : Wood & Chan notation; 
  #since we are going upto n and not stopping at (n-1), the condition is modified!
  return("g" = g)
}

# forming the circulant matrix:
circulant=function(x){
  n = length(x)
  mat = matrix(0, n, n)
  for (j in 1:n) {
    mat[j, ] <- c(x[-(1:(n+1-j))], x[1:(n+1-j)])
  }
  return(mat)
}

# Function for forming the vector of circulant matrix:
circ_vec=function(knot,g,nu,l,tausq){
  delta_N=1/(length(knot)-1)
  m=2**g
  cj=integer()
  for(j in 1:m){
    if(j<=(m/2))
      cj[j]=(j-1)*delta_N
    else
      cj[j]=(m-(j-1))*delta_N
  }
  x=(tausq*MK(cj,0,l,nu))
  return(x)
}

# Function for finding a g such that C is nnd:
eig.eval=function(knot,g,nu,l,tausq){
  vec=circ_vec(knot,g,nu,l,tausq)
  C=circulant(vec)
  ev=min(eigen(C)$values)
  return(list("vec" = vec, "min.eig.val" = ev))
}

# Function for finding a g such that C is nnd:
# without forming the circulant matrix and without computing eigen values:
C.eval=function(knot,g,nu,l,tausq){
  vec=circ_vec(knot,g,nu,l,tausq)
  val=fft(vec) # eigenvalues will be real as the circulant matrix formed by the 
  # vector is by construction is symmetric!
  ev=min(Re(val))
  return(list("vec" = vec, "min.eig.val" = ev))
}


nnd_C=function(knot,g,nu,l,tausq){
  C.vec=C.eval(knot,g,nu,l,tausq)$vec
  eval=C.eval(knot,g,nu,l,tausq)$min.eig.val
  if(eval>0)
    return(list("cj" = C.vec,"g" = g))
  else{
    g=g+1
    nnd_C(knot,g,nu,l,tausq)
  }
}

# computing the eigen values of C using FFT:
eigval=function(knot,nu,l,tausq){
  g=min_g(knot)
  c.j=nnd_C(knot,g,nu,l,tausq)$cj
  lambda=Re(fft(c.j))
  if(min(lambda)>0)
    return(lambda)
  else
    stop("nnd condition is NOT satisfied!!!")
}


#################################################################
########## Samples drawn using Wood and Chan Algorithm ##########
#################################################################
samp.WC=function(knot,nu,l,tausq){
  N=length(knot)
  lambda=eigval(knot,nu,l,tausq)
  m=length(lambda)
  samp.vec=rep(0,N)
  a=rep(0,m)
  a[1]=sqrt(lambda[1])*rnorm(1)/sqrt(m)
  a[(m/2)+1]=sqrt(lambda[(m/2)+1])*rnorm(1)/sqrt(m)
  i=sqrt(as.complex(-1))
  for(j in 2:(m/2)){
    uj=rnorm(1); vj=rnorm(1)
    a[j]=(sqrt(lambda[j])*(uj + i*vj))/(sqrt(2*m))
    a[m+2-j]=(sqrt(lambda[j])*(uj - i*vj))/(sqrt(2*m))
  }
  samp=fft(a)
  samp.vec=Re(samp[1:N])
  return(samp.vec)
}

#############################################
########## Functions for using ESS ##########
#############################################
ESS = function(beta,nu_ess,y,X,sigsq,cj,x1){
  thetamin = 0; 
  thetamax = 2*pi;
  
  theta = runif(1,thetamin,thetamax); 
  thetamin = theta - 2*pi; 
  thetamax = theta;
  betaprime = beta*cos(theta) + nu_ess*sin(theta)
#  alpha=min(1,(lik(y,X,sigsq,betaprime,cj,x1)/lik(y,X,sigsq,beta,cj,x1)))
 # alpha=min(1,exp(loglik(y,X,sigsq,betaprime)*((min(betaprime) >= 0) && ((x1 + sum(cj*betaprime))<= 0))-loglik(y,X,sigsq,beta)))
  if( (min(betaprime) >= 0 && (x1 + sum(cj*betaprime))<= 0)){
    alpha=min(1,exp(loglik(y,X,sigsq,betaprime)-loglik(y,X,sigsq,beta)))
  }else{
    alpha = 0
  }
  u = runif(1)
  
  if(alpha > u){
    return(betaprime)
  }else{
    while(alpha < u){
      if(theta < 0)
        thetamin = theta
      else
        thetamax = theta
      theta = runif(1,thetamin,thetamax)
      betaprime = beta*cos(theta) + nu_ess*sin(theta)
     # alpha=min(1,(lik(y,X,sigsq,betaprime,cj,x1)/lik(y,X,sigsq,beta,cj,x1)))
      if(min(betaprime) >= 0  &&  x1 + sum(cj*betaprime) <= 0){
        alpha=min(1,exp(loglik(y,X,sigsq,betaprime)-loglik(y,X,sigsq,beta)))
      }else{
        alpha = 0
      }
      
    }
    return(betaprime)
  }
}

######################################################
########## Functions for using ESS (proton) ##########
######################################################
ESS_proton = function(beta,nu_ess,ytilde,X,sigsq,tausq,Inv_tr,Inv_cov_br_2,cj,x1){
  thetamin = 0; 
  thetamax = 2*pi;
  
  theta = runif(1,thetamin,thetamax); 
  thetamin = theta - 2*pi; 
  thetamax = theta;
  betaprime = beta*cos(theta) + nu_ess*sin(theta)
  alpha=min(1,(lik_beta(ytilde,X,sigsq,tausq,Inv_tr,Inv_cov_br_2,betaprime,cj,x1)/lik_beta(ytilde,X,sigsq,tausq,Inv_tr,Inv_cov_br_2,beta,cj,x1)))
  u = runif(1)
  
  if(alpha > u){
    return(betaprime)
  }else{
    while(alpha < u){
      if(theta < 0)
        thetamin = theta
      else
        thetamax = theta
      theta = runif(1,thetamin,thetamax)
      betaprime = beta*cos(theta) + nu_ess*sin(theta)
      alpha=min(1,(lik_beta(ytilde,X,sigsq,tausq,Inv_tr,Inv_cov_br_2,betaprime,cj,x1)/lik_beta(ytilde,X,sigsq,tausq,Inv_tr,Inv_cov_br_2,beta,cj,x1)))
    }
    return(betaprime)
  }
}

################################################################

ESS2_proton = function(beta,nu_ess,ytilde,X,sigsq,tausq,Inv_tr,Inv_cov_br_2,cj,x1){
  
  thetamin = 0; 
  thetamax = 2*pi;
  
  theta = runif(1,thetamin,thetamax); 
  thetamin = theta - 2*pi; 
  thetamax = theta;
  betaprime = beta*cos(theta) + nu_ess*sin(theta)
  if (min(betaprime) >= 0  &&  x1 + sum(cj*betaprime) <= 0){
    alpha=min(1,exp(loglik_beta(ytilde,X,sigsq,tausq,Inv_tr,Inv_cov_br_2,betaprime,x1)-loglik_beta(ytilde,X,sigsq,tausq,Inv_tr,Inv_cov_br_2,beta,x1)))
  }else{alpha = 0}
  u = runif(1)
  i = 0 
  if(alpha > u){
    betanew = (betaprime)
  }else{
    
    while(alpha < u){
      if(theta < 0)
        thetamin = theta
      else
        thetamax = theta
      theta = runif(1,thetamin,thetamax)
      betaprime = beta*cos(theta) + nu_ess*sin(theta)
      #alpha=min(1,(lik_beta(ytilde,X,sigsq,tausq,Inv_tr,Inv_cov_br_2,betaprime,cj,x1)/lik_beta(ytilde,X,sigsq,tausq,Inv_tr,Inv_cov_br_2,beta,cj,x1)))
      if (min(betaprime) >= 0  &&  x1 + sum(cj*betaprime) <= 0){
        alpha=min(1,exp(loglik_beta(ytilde,X,sigsq,tausq,Inv_tr,Inv_cov_br_2,betaprime,x1)-loglik_beta(ytilde,X,sigsq,tausq,Inv_tr,Inv_cov_br_2,beta,x1)))
      }else{alpha = 0}
      
     # i = i+1
     # print(i)
    }
    betanew = (betaprime)
  }
  
  
}

#################################################


ESS_beta0 = function(x1,nu_x1,m,s,cj,xi){
  thetamin = 0; 
  thetamax = 2*pi;
  
  theta = runif(1,thetamin,thetamax); 
  thetamin = theta - 2*pi; 
  thetamax = theta;
  x1prime = x1*cos(theta) + nu_x1*sin(theta)
  if (x1prime + sum(cj*xi) <= 0){
  alpha=min(1,exp(loglik_beta0(x1prime,m,s)-loglik_beta0(x1,m,s)))
  }else{alpha = 0}
  
 # alpha=min(1,(lik_beta0(x1prime,m,s,cj,xi)/lik_beta0(x1,m,s,cj,xi)))
 
  u = runif(1)
  # u = log(runif(1))
  if(alpha > u){
    return(x1prime)
  }else{
    while(alpha < u){
      if(theta < 0)
        thetamin = theta
      else
        thetamax = theta
      theta = runif(1,thetamin,thetamax)
      x1prime = x1*cos(theta) + nu_x1*sin(theta)
      
     # alpha=min(1,(lik_beta0(x1prime,m,s,cj,xi)/lik_beta0(x1,m,s,cj,xi)))
      
      if (x1prime + sum(cj*xi) <= 0){
        alpha=min(1,exp(loglik_beta0(x1prime,m,s)-loglik_beta0(x1,m,s)))
      }else{alpha = 0}
      
    }
    return(x1prime)
  }
}

ESS_xi0 = function(xi0,nu_xi0,m0,s0,sig0){
  thetamin = 0; 
  thetamax = 2*pi;
  
  theta = runif(1,thetamin,thetamax); 
  thetamin = theta - 2*pi; 
  thetamax = theta;
  xi0prime = xi0*cos(theta) + nu_xi0*sin(theta)
  alpha=min(1,(lik_xi0(xi0prime,m0,s0,sig0)/lik_xi0(xi0,m0,s0,sig0)))
  # alpha=min(0,(loglik_beta0(x1prime,m,s,cj,xi)-loglik_beta0(x1,m,s,cj,xi))*(((x1 + sum(cj*xi))<=0) && ((x1prime + sum(cj*xi))<=0)))
  u = runif(1)
  # u = log(runif(1))
  if(alpha > u){
    return(xi0prime)
  }else{
    while(alpha < u){
      if(theta < 0)
        thetamin = theta
      else
        thetamax = theta
      theta = runif(1,thetamin,thetamax)
      xi0prime = xi0*cos(theta) + nu_xi0*sin(theta)
      alpha=min(1,(lik_xi0(xi0prime,m0,s0,sig0)/lik_xi0(xi0,m0,s0,sig0)))
      # alpha=min(0,(loglik_beta0(x1prime,m,s,cj,xi)-loglik_beta0(x1,m,s,cj,xi))*((x1 + sum(cj*xi))<=0 && (x1prime + sum(cj*xi))<=0))
    }
    return(xi0prime)
  }
}


## Defining the loglik function to be used in ESS:
## loglik calculates the log of the likelihood:
lik=function(y,X,sigsq,beta,cj,x1){
  mu=y-(X%*%beta)
  val=exp(-sum(mu^2)/(2*sigsq))*((min(beta) >= 0) && ((x1 + sum(cj*beta))<= 0))
  return(val)
}

loglik=function(y,X,sigsq,beta){
  mu=y-(X%*%beta)
  val=-sum(mu^2)/(2*sigsq)
  return(val)
}

# Inverse of covariance matrix 
#Invmat = function(Covmatp){
#  Sigma = Covmatp[-1,-1]-Covmatp[-1,1]%*%Covmatp[1,-1]/Covmatp[1,1]
#  Inv_cov = solve(Sigma)
#  Inv_covsub = solve(Covmatp[-1,-1])
#  return(list(Sigma, Inv_cov, Inv_covsub))
#}

## Defining the loglik function to be used in ESS for joint \xi:
## loglik calculates the log of the likelihood:
lik_beta=function(ytilde,X,sigsq,tausq,Inv_tr,Inv_cov_br_2,beta,cj,x1){
  # Inv_tr:top-right block of the inverse of covariance matrix
  mu=ytilde-(X%*%beta)
 # mu1=Covmatp[-1,1]*x1/Covmatp[1,1]
  #val=exp(-sum(mu^2)/(2*sigsq)-x1*sum(Inv_tr*beta)/tausq+beta%*%Inv_cov_br_2%*%beta/2/tausq)*((min(beta) >= 0) && ((x1 + sum(cj*beta))<= 0))
  val=exp(-sum(mu^2)/(2*sigsq)-x1*sum(as.numeric(Inv_tr)*beta)/tausq+beta%*%Inv_cov_br_2%*%beta/2/tausq)*((min(beta) >= 0) && ((x1 + sum(cj*beta))<= 0))
  return(val)
}

loglik_beta=function(ytilde,X,sigsq,tausq,Inv_tr,Inv_cov_br_2,beta,x1){
  # Inv_tr:top-right block of the inverse of covariance matrix
  mu=ytilde-(X%*%beta)
  # mu1=Covmatp[-1,1]*x1/Covmatp[1,1]
  #val=exp(-sum(mu^2)/(2*sigsq)-x1*sum(Inv_tr*beta)/tausq+beta%*%Inv_cov_br_2%*%beta/2/tausq)*((min(beta) >= 0) && ((x1 + sum(cj*beta))<= 0))
  val=-sum(mu^2)/(2*sigsq)-x1*sum(as.numeric(Inv_tr)*beta)/tausq+beta%*%Inv_cov_br_2%*%beta/2/tausq
  return(val)
}

loglik_beta0=function(x1,m,s){
#  val=m*x1/s^2
  val=-(m-x1)^2/s^2/2+x1^2/2
  return(val)
}

lik_beta0=function(x1,m,s,cj,xi){
  val=exp(-(m-x1)^2/s^2/2+x1^2/2)*((x1 + sum(cj*xi))<=0)
  return(val)
}

lik_xi0=function(xi0,m0,s0,sig0){
  val=exp(-(xi0-m0)^2/s0^2/2+xi0^2/2)*(xi0>=1-sig0 && xi0<=1+sig0)
}








#######################
### Model Checking ####
#######################

#### WAIC ####
getWAIC = function(y,fhat_sam,fmean,sig_sam){
  
  n = length(y)
  n_sam = length(sig_sam)
  
  # likelihood function 
  #mu_sam = xi0_sam + xi1_sam*x + X%*%xi_sam
  lik.fun = function(y,mu,sig) dnorm(y, mean = mu, sd = sqrt(sig), log = F)
  lik.value = matrix(nrow = n, ncol = n_sam)
  
  # likelihood values
  for(i in 1:n){
    for(j in 1:n_sam){
      lik.value[i,j] = lik.fun(y[i], fhat_sam[i,j], sig_sam[j])
    }
  }
  
  # lik.value the likelihood function evaluated over each observation y_i and mcmc samples theta_s.
  lppd = sum(log(rowMeans(lik.value)))
  
  # P_waic_1
  P_waic1 = 2*sum(log(rowMeans(lik.value)) - rowMeans(log(lik.value)))
  
  # P_waic_2 
  lik.avg = lik.temp = c()
  mu.avg = mean(fmean)
  sig.avg = mean(sig_sam)
 # for (i in 1:n){lik.avg[i] = lik.fun(y[i], mu.avg, sig.avg)}
  for (i in 1:n){lik.temp[i] = sum((log(lik.value[i,])-mean(log(lik.value[i,])))^2)/(n_sam-1)}
  P_waic2 = sum(lik.temp)
  
  # ellppd
  elppd_WAIC1 = -2*(lppd - P_waic1)
  elppd_WAIC2 = -2*(lppd - P_waic2)
  
  return(list(elppd_WAIC1, elppd_WAIC2))
  
}


#### DIC ####
getDIC = function(y,fhat_sam,fmean,sig_sam){
  pDIC = c()
  n = length(y)
  n_sam = length(sig_sam)
  
  # likelihood function 
  #mu_sam = xi0_sam + xi1_sam*x + X%*%xi_sam
  lik.fun = function(y,mu,sig) dnorm(y, mean = mu, sd = sqrt(sig), log = F)
  lik.value = matrix(nrow = n, ncol = n_sam)
  # likelihood values
  for(i in 1:n){
    for(j in 1:n_sam){
      lik.value[i,j] = lik.fun(y[i], fhat_sam[i,j], sig_sam[j])
    }
  }
  sig.avg = mean(sig_sam) 
  lik.avg=c()
  for(i in 1:n){
    lik.avg[i] = lik.fun(y[i],fmean[i],sig.avg)
    pDIC[i] = 2*(log(lik.avg[i])-mean(lik.value[i,]))
  }
  
  DIC = -2*log(lik.avg)+2*pDIC
  
  return(DIC)
}


##### BIC #####
getBIC = function(y,fhat_sam,fmean,xi_sam,sig_sam){
  n = length(y)
  n_sam = length(sig_sam)
  N = dim(xi_sam)[2]
  k = N+5
  # likelihood function 
  #mu_sam = xi0_sam + xi1_sam*x + X%*%xi_sam
  lik.fun = function(y,mu,sig) dnorm(y, mean = mu, sd = sqrt(sig), log = F)
  lik.value = matrix(nrow = n, ncol = n_sam)
  sig.avg = mean(sig_sam) 
  fmean.avg = mean(fmean) 
  lik.avg=BIC=c()
  for(i in 1:n){
    lik.avg[i] = lik.fun(y[i],fmean.avg,sig.avg)
    #-2*log(p(y|\hat{\theta})+k\log(n)  
    BIC[i] = -2*log(lik.avg[i])+k*log(n)
  }
  return(BIC)
}



