# wrappers to unify probability estimation interface and output ----------------

sov = function(mu, Sigma, lb, ub, params) {
  gb <- tlrmvnmvt::GenzBretz(params$N)
  
  prob <- tlrmvnmvt::pmvn(lb, ub, mu, Sigma, algorithm = gb)
  result <- data.frame(variable = c("estimate", "error"), 
                       value = c(prob, attr(prob, "error")))
  return(result)
}



