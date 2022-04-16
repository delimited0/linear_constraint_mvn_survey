#####################################################################
########## Function for MCMC samples using ESS and WC algo ##########
#####################################################################

library(MASS); library(FastGP)

### Function for drawing posterior samples using ESS with fixed hyperparameters:

########### For convex functions estimation ##############
RobustGP=function(Y,X,sig,nu,N,l,sig0,prob,mcmc,brn,thin,tau.in,sig.in,xi0.in,xi1.in,xi.in,
                 xi0.fix,xi1.fix,tau.fix,sig.fix,xi.fix,sseed,
                 verbose,return.plot,return.traplot){
  # y: Response variable (G_E) 
  # x: Covariates (Q^2)
  # sig: noise level 
  # nu: smoothness parameter of Matern
  # N: the number of knots - 1
  # l: length-scale parameter of Matern
  # prob: a scalar [0, 1] specifying the probability within the credible interval
  # mcmc, brn, thin : mcmc samples, burning and thinning for MCMC
  # tau.in, sig.in, xi0.in, xi1.in, xi.in : initial values (supplied by user or use the default value)
  # xi0.fix, xi1.fix, tau.fix, sig.fix, xi.fix : if fixed values are to use
  # verbose : logical; if TRUE, prints currnet status; default is TRUE
  # return.plot : logical; if true a plot of estimate with 95% CI is returned; default is TRUE
  # return.traplot : logical; if true traceplots are returned; default is TRUE
  
  # OUTPUT: 
  # (1) r_est: estimate of proton redius;  
  # (2) r_CI: 95% credible interval 
  # (2) data: data matrix with x, ytrue, y 
  # (3) post.sam: Posterior samples on xi,xi0,xi1,tau,sig and fhat with posterior mean, 95% CI of fhat
  # (4) model.check: a list object including WAIC; DIC; BIC
  # (5) HD_rp: Highest density set of proton redius
  # (6) r_Ess: Effective number of parameters for mcmc 
  # (7) r_mcse: the standard error for mcmc samples of proton redius 
  
  Y = as.numeric(as.matrix(Y))
  X = as.numeric(as.matrix(X))
  
  n = length(X)
  xmax = max(X)
  x = X/xmax
  if(!missing(sig)){
  y = Y + sig*rnorm(n,0,1)
  }else{y = Y}
  
  post.sam = ESS.joint(y, x, nu=nu, N=N, l=l,
                       sig0=sig0,
                       sseed=sseed, mcmc=mcmc, brn=brn, thin=thin,
                       xi0.fix = xi0.fix,
                       return.plot=return.plot,return.traplot=return.traplot)
  
  # (1) Estimation and fitted curve:
  fmean=post.sam$fmean
  f_low=post.sam$f_low
  f_upp=post.sam$f_upp
  
  ub=max(Y,f_low,f_upp,fmean)
  lb=min(Y,f_low,f_upp,fmean)

 # if(return.plot){  
#  df = data.frame(x,fmean)
#  G <- ggplot(df, aes(x=x, y=fmean), color=variable)+ 
#    geom_line(aes(x=x, y=fmean), colour="blue")+
#    geom_point(aes(x=x, y=Y), colour="red",cex = 0.3)+
#    geom_ribbon(aes(ymin=f_low, ymax=f_upp), alpha=0.8)
  
#  G <- G + theme_bw()
#  Glabs<- G+labs(x = "Q^2", y = "G_E")
#  Glabs + theme( 
#    axis.title.x = element_text(size=8, face="bold"),
#    axis.title.y = element_text(size=8, face="bold"),
#    legend.text = element_text(colour="blue", size=10, face="bold"))
#  if(return.zoominplot){
#  Gz <- G + coord_cartesian(xlim = c(0,0.005), ylim = c(0.998,1.002), expand = TRUE)
#  Gzlabs<- Gz+labs(x = "Q^2", y = "G_E")
#  Gzlabs2<- Gzlabs + theme( 
#    axis.title.x = element_text(size=8, face="bold"),
#    axis.title.y = element_text(size=8, face="bold"),
#    legend.text = element_text(colour="blue", size=10, face="bold"))
#  vp <- viewport(width = 0.45, height = 0.45, x = 0.75, y = 0.75)
#  print(Gzlabs2, vp=vp)
#  }
  
  
#}
  if(return.plot){
    library(ggplot2)
    df = data.frame(x,fmean)
    G <- ggplot(df, aes(x=x, y=fmean), color=variable)+ 
      geom_line(aes(x=x, y=fmean), colour="blue")+
      geom_point(aes(x=x, y=Y), colour="red",cex = 0.3)+
      geom_ribbon(aes(ymin=f_low, ymax=f_upp), alpha=0.8)
    G
    G <- G + theme_bw()
    Glabs<- G+labs(x = "Q^2", y = "G_E")
    Glabs + theme( 
      axis.title.x = element_text(size=8, face="bold"),
      axis.title.y = element_text(size=8, face="bold"),
      legend.text = element_text(colour="blue", size=10, face="bold"))
  }
  
  
  
  ### (2) model checking ###
  ### mcmc samples ##
  fhat_sam=post.sam$fhat_sam
  fmean=post.sam$fmean 
  xi_sam=post.sam$xi_sam
  xi0_sam=post.sam$xi0_sam
  xi1_sam=post.sam$xi1_sam
  sig_sam=post.sam$sig_sam
  tau_sam=post.sam$tau_sam
  
  WAIC = getWAIC(y,fhat_sam,fmean,sig_sam)         # joint
  DIC = getDIC(y,fhat_sam,fmean,sig_sam)           # pointwise--vector
  BIC = getBIC(y,fhat_sam,fmean,xi_sam,sig_sam)    # pointwise--vector
  
  model.check = list(WAIC,DIC,BIC)
  
  ### (3) Highest Density Interval
  # For proton radius and normalizing factor
  r_sam = sqrt(-6*xi1_sam/xi0_sam/xmax)
 # HDprob = ifelse(missing(HDprob),0.8, HDprob)
  if(missing(prob)) {prob=0.95}
  HD_rp = hdi(r_sam, credMass=prob)
  HD_xi0 = hdi(xi0_sam, credMass=prob)
  

  ### (4) standard error for mcmc samples ####
  ## effiective sample size for mcmc sample of the proton radius 
  
  r_sam = sqrt(-6*xi1_sam/xi0_sam/xmax)
  r_Ess = mcmcse::ess(r_sam)
  r_mcse = mcmcse::mcse(r_sam,method="obm")
  
  
  ### (5) Final estimation !! ###
  r_est = median(sqrt(-6*xi1_sam/xi0_sam/xmax)) 
  r_CI = quantile(sqrt(-6*xi1_sam/xi0_sam/xmax), probs=c(0.025,0.975)) 
  
  data = data.frame(x,y)
  
  Results = list(r_est,r_CI,r_sam,data,post.sam,model.check,HD_rp,r_Ess,r_mcse)


  return(Results)
  
  }


#END