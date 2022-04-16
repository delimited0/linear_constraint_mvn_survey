#### GP covariance kernel -- matern kernel

maternCov = function(u, l, nu){
  M = function(x, y ,l, nu){
    ifelse(abs(x-y)>0, (sqrt(2*nu)*abs(x-y)/l)^nu/(2^(nu-1)*gamma(nu))*besselK(x=abs(x-y)*sqrt(2*nu)/l, nu=nu), 1.0)
  } 
  
  M = function(x, y ,l, nu){
    ifelse(abs(x-y)>0, (sqrt(2*nu)*abs(x-y)/l)^nu/(2^(nu-1)*gamma(nu))*besselK(x=abs(x-y)*sqrt(2*nu)/l, nu=nu), 1.0)
  } 
  
  M_x = function(x, y, l, nu) -2^(3/2-nu/2)/gamma(nu)*nu*(x-y)/(l^2)*(sqrt(nu)/l*abs(x-y))^(nu-1) *besselK(x=abs(x-y)*sqrt(2*nu)/l, nu=nu-1) 
  
  M_xx = function(x, y, l, nu) { 2^(1/2-nu/2)*(nu/(l^2))^(3/2)/gamma(nu)*(sqrt(nu)/l*abs(x-y))^(nu-2)*
      (2*sqrt(nu)/l*(x-y)^2*besselK(x=abs(x-y)*sqrt(2*nu)/l, nu=nu-2)-sqrt(2)*abs(x-y)*besselK(x=abs(x-y)*sqrt(2*nu)/l, nu=nu-1))}
  
  M_xy = function(x, y, l, nu) - M_xx(x, y, l, nu)
  
  M_xxy = function(x, y, l, nu) {2^(-nu/2)*nu/l^2*(x-y)*(sqrt(nu)/l*abs(x-y))^nu/(gamma(nu)*abs(x-y)^3)* 
      (4*sqrt(2)*sqrt(nu)/l*(x-y)^2*besselK(x=abs(x-y)*sqrt(2*nu)/l, nu=nu-3) -
         12*abs(x-y)*besselK(x=abs(x-y)*sqrt(2*nu)/l, nu=nu-2))}
  
  
  M_xxyy = function(x, y, l, nu) {2^(2-nu/2)*nu*(sqrt(nu)/l*abs(x-y))^(nu-1)/(l^5*gamma(nu)*abs(x-y)^3)*
      (3*sqrt(2)*l*abs(x-y)*(l^2*(nu-3)-2*nu*(x-y)^2)*besselK(x=abs(x-y)*sqrt(2*nu)/l, nu=nu-3)+
         sqrt(nu)*(3*l^2*(x-y)^2+2*nu*(x-y)^4)*besselK(x=abs(x-y)*sqrt(2*nu)/l, nu=nu-4))}
  
  
  N = length(u)
  mcov = matrix(nrow = N+2, ncol = N+2)
  mcov_off = c()
  mcov[1,1] = M(0,0,l,nu)
  mcov[2,1] = 0
  mcov[1,2] = mcov[2,1]
  mcov[2,2] = M_xy(10^(-8),0,l,nu)
  mcov[3,1] = M_xx(10^(-8),0,l,nu)
  mcov[1,3] = mcov[3,1]
  mcov[3,2] = M_xxy(10^(-8),0,l,nu)
  mcov[2,3] = mcov[3,2]
  mcov[-c(1:3),1] = M_xx(u[-1], 0, l, nu) 
  mcov[1,-c(1:3)] = t(mcov[-c(1:3),1])
  mcov[-c(1:3),2]=  M_xxy(u[-1], 0, l, nu) 
  mcov[2,-c(1:3)] = t(mcov[-c(1:3),2])
  mcov_off0 = 12*nu^4/((2-3*nu+2*nu^2)*l^4)
  for(k in 1:(N-1)){ mcov_off[k] = M_xxyy(k*dN, 0, l, nu) }
  
  for( i in 1:length(u)){
    
    mcov[i+2,i+2] = mcov_off0 
    
    for(j in 1:length(u)){
      
      for(k in 1:(N-1)){
        if(abs(i-j) == k)  mcov[i+2,j+2] = mcov_off[k]
      }
    }
  }
  
  return(mcov)
  
}

