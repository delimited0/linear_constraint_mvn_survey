#####################################################################
########## Function for MCMC samples using ESS and WC algo ##########
#####################################################################

library(MASS); library(FastGP)

### Function for drawing posterior samples using ESS with fixed hyperparameters:

########### For convex functions estimation ##############
ESS.joint = function(y, x, nu, N, l, sig0, 
                     mcmc, brn, thin, 
                     tau.in, sig.in, xi0.in, xi1.in, xi.in,
                     xi0.fix, xi1.fix, tau.fix, sig.fix, xi.fix,
                     sseed, verbose, return.plot, return.traplot) {
  # y:Response variable; x: vector to form design matrix \Psi (n X N+1)
  # N: (N+1) is the number of knots
  # nu:smoothness parameter of Matern; l:length-scale parameter of Matern
  # nt & a: the parameters for numerical calculation of the integration of matern kernel
  # mcmc, brn, thin : mcmc samples, burning and thinning for MCMC
  # tau.in, sig.in, xi0.in, xi1.in, xi.in : initial values (supplied by user or use the default value)
  # xi0.fix, xi1.fix, tau.fix, sig.fix, xi.fix : if fixed values are to use
  # verbose : logical; if TRUE, prints currnet status; default is TRUE
  # return.plot : logical; if true a plot of estimate with 95% CI is returned; default is TRUE
  # return.traplot : logical; if true traceplots are returned; default is TRUE
  
  # OUTPUT: Posterior sampless on xi,xi0,xi1,tau,sig and fhat with posterior mean, 95% CI of fhat
  ptm = proc.time()
  
  if (length(y) != length(x))
    stop("y and x should be of same length!")
  n = length(y)
  if (missing(N))
    N = ceiling(n/2) - 1
  int_length = 1/N
  my_knots = seq(0, 1, by=int_length)
  X = des.mat2(x, my_knots, int_length)
  cj = psi_j(1, my_knots, int_length)
  
  if (missing(nu))
    stop("nu needs to be supplied")
  if (nu==0)
    stop("nu cannot be zero")
  if (!missing(l) && l==0)
    stop("l cannot be zero")
  
  if (missing(return.plot))
    return.plot = TRUE
  if (missing(return.traplot))
    return.traplot = TRUE
  
  if(!missing(sseed))
    set.seed(sseed)
  if(missing(sseed))
    set.seed(Sys.Date())
  
  if(missing(verbose))
    verbose=TRUE
  
  if(missing(l)) 
    l=l_est(nu,c(my_knots[1],my_knots[length(my_knots)]),0.05)
  
  # prior covariance K:
  K = covmat(my_knots,nu,l)
  # Joint covariance:
  nt = 1000  # parameters for numerical cal 
  Covmatp = covmat_joint(my_knots, l, nu, nt, a=20)
  # precision:
  Inv_cov = solve(Covmatp+10^(-8)*diag(dim(Covmatp)[1])) # inverse of whole
  # blocked matrices 
  dr = dim(Covmatp)[1]
  Cov_tl = Covmatp[1,1] #C11
  Cov_br = Covmatp[c(2:dr),c(2:dr)] #C22
  Cov_bl = Covmatp[c(1:(dr-1)),1] #C21
  Cov_tr = t(Cov_bl) #C12
  Inv_tl_1 = 1/Cov_tl #S11_1
  Inv_br = tinv(Cov_br-Cov_bl%*%Cov_tr*as.numeric(Inv_tl_1)) #S22
  Inv_tr = -Cov_bl%*%Inv_br*Inv_tl_1 #S12
  Inv_tl_2 = Inv_tl_1*Cov_tr%*%Inv_br%*%Cov_bl*Inv_tl_1 #S11_2
  Inv_tl = Inv_tl_1 + Inv_tl_2 #S11
  Inv_bl = t(Inv_tr) #S21
  # to use the WC algorithm
  Inv_cov_br_1 = tinv(Cov_br)
  Inv_cov_br_2 = -Inv_cov_br_1%*%Cov_bl%*%Cov_tr%*%Inv_cov_br_1/as.numeric(-Inv_tl_1+Cov_tr%*%Inv_cov_br_1%*%Cov_bl)
  
  if(missing(mcmc))
    mcmc=5000
  if(missing(brn))
    brn=1000
  if(missing(thin))
    thin=1
  em=mcmc+brn
  ef=mcmc/thin
  
  if(!missing(tau.fix))
    tau.in=tau.fix
  if(!missing(sig.fix))
    sig.in=sig.fix
  if(!missing(xi0.fix))
    xi0.in=xi0.fix
  if(!missing(xi1.fix))
    xi1.in=xi1.fix
  if(!missing(xi.fix))
    xi.in=xi.fix
  
  if(missing(tau.fix) && missing(tau.in))
    tau.in=1
  if(missing(sig.fix) && missing(sig.in))
    sig.in=1
  if(missing(xi0.fix) && missing(xi0.in))
    xi0.in=1
  if(missing(xi.fix) && missing(xi.in))
  # selecting the initial values 
  post.sam1 =  ESS.pro(y,x,nu=nu,N=N,l=l,sseed=67,mcmc=2000,brn=500,thin=5,
          return.plot = FALSE,return.traplot = FALSE)
  #xi1.in=median(post.sam1$xi1_sam)
  xi.in=rowMeans(post.sam1$xi_sam)
  if(missing(xi1.fix) && missing(xi1.in))
    xi1.in=-sum(cj*xi.in)-0.1
    
  med_x1 = "ESS"
  med_xi0 = "ESS"
  
  tau=tau.in
  sig=sig.in
  xi0=xi0.in
  xi1=xi1.in
  xi_in=xi.in
  xi1_in=xi1.in
  xi0_in=xi0.in
  
  
  xi_sam=matrix(0,N+1,ef)
  xi0_sam=rep(0,ef)
  xi1_sam=rep(0,ef)
  tau_sam=rep(0,ef)
  sig_sam=rep(0,ef)
  fhat_sam=matrix(0,n,ef)
  
  if(verbose)
    print("MCMC sample draws:")
  
  for(i in 1:em){
    # sampling \Xi:
    if(missing(xi.fix)){
    y_tilde = y - xi0_in - xi1_in*x
    nu.ess = as.vector(samp.WC(my_knots,nu,l,tau))
    xi_out = ESS2_proton(xi_in,nu.ess,y_tilde,X,sig,tau,Inv_tr,Inv_cov_br_2,cj,xi1)  #(beta,nu_ess,y,X,sigsq,cj,x1)
    #xi_out = pmax(0,xi_out)
}else{
      xi_out = xi_in
    }
    # sampling \xi_0:
    Xxi = as.vector(X %*% xi_out)
    y_star = y - xi1*x - Xxi
    m0 = mean(y_star)
    s0 = sqrt(sig/n) 
    if(missing(xi0.fix)){
      if(med_xi0 == "Inv"){
         lb = pnorm(1-sig0,m0,s0)
         ub = pnorm(1+sig0,m0,s0)
         u0 = runif(1,lb,ub)
         xi0 = m0+s0*qnorm(u0)
  }else if(med_xi0 == "ESS"){
         nu_xi0 = rnorm(1,0,1)
         xi0 = ESS_xi0(xi0_in,nu_xi0,m0,s0,sig0)
  }
    }else{xi0 = xi0.fix}    
      
    
    # sampling \xi_1:
    y1 = y - xi0 - Xxi
    if(missing(xi1.fix)){ #draw from univariate truncated normal
    s = 1/sqrt(sum(x^2)/sig+as.numeric(Inv_tl)/tau)
    m = (sum(x*y1)/sig-sum(Inv_tr*xi_out)/tau)*s^2
    #ESS 
    if(med_x1 == "ESS"){
    nu_x1 = rnorm(1,0,1)
    xi1=ESS_beta0(xi1_in,nu_x1,m,s,cj,xi_out)
    }else if(med_x1 == "Inv"){
    # Inversion
    u=runif(1,0,pnorm(-sum(cj*xi_out),m,s))
    xi1=m+s*qnorm(u)
    }
    }
    
    # sampling \sigma^2:
    y0 = y - xi0 - xi1*x - Xxi
    if(missing(sig.fix))
      sig = 1/rgamma(1,shape = n/2,rate = sum(y0^2)/2)
    
    # sampling \tau^2:
    if(missing(tau.fix))
      xi_joint = c(xi1,xi_out)
      tau = 1/rgamma(1, shape = (N+2)/2, rate = (t(xi_joint)%*%Inv_cov%*%xi_joint)/2)
    
    # storing MCMC samples:
    if(i > brn && i%%thin == 0){
      browser()
      xi_sam[,(i-brn)/thin]=xi_out
      xi0_sam[(i-brn)/thin]=xi0
      xi1_sam[(i-brn)/thin]=xi1
      sig_sam[(i-brn)/thin]=sig
      tau_sam[(i-brn)/thin]=tau
      fhat_sam[,(i-brn)/thin]=xi0 + xi1*x + Xxi
    }
    
    if(i%%2000==0 && verbose){
      print(i)
    }
    
    # renewing the initial value:
    xi_in = xi_out
    xi1_in = xi1
    xi0_in = xi0
  }
  
  if(return.traplot){
    library(mcmcplots)
    mx=min(N+1,25)
    xi_mat=matrix(xi_sam[1:mx,],nrow = mx,
                  dimnames = list(c(paste("xi[",1:mx,"]",sep = "")),NULL))
    traplot(t(xi_mat))
    par(mfrow=c(1,1))
    traplot(as.mcmc(xi0_sam),main=expression(paste(xi[0])))
    traplot(as.mcmc(xi1_sam),main=expression(paste(xi[1])))
    traplot(as.mcmc(sig_sam),main=expression(paste(sigma^2)))
    traplot(as.mcmc(tau_sam),main=expression(paste(tau^2)))
  }
  
  fmean=rowMeans(fhat_sam)
  q=apply(fhat_sam,1,function(x) return(quantile(x,c(0.025,0.975))))
  f_low=as.vector(q[1,])
  f_upp=as.vector(q[2,])
  ub=max(f_low,f_upp,fmean)
  lb=min(f_low,f_upp,fmean)
  
  if(return.plot){
    library(ggplot2)
    par(mfrow=c(1,1))
    plot(x,fmean,type="l",lwd=2,main="Black: Estimate, Blue: 95% CI, Red: Data",ylim=c(lb,ub),ylab = "G_E",xlab = "Q^2")
    points(x,y,cex=0.5,col="red")
    lines(x,f_low,lwd=2,lty=2,col="blue")
    lines(x,f_upp,lwd=2,lty=2,col="blue")
  }
  
  tm=proc.time()-ptm
  return(list("time"=tm,"xi_sam"=xi_sam,"xi0_sam"=xi0_sam,"xi1_sam"=xi1_sam,"sig_sam"=sig_sam,
              "tau_sam"=tau_sam,"fhat_sam"=fhat_sam,"fmean"=fmean,"f_low"=f_low,"f_upp"=f_upp))
}


#END