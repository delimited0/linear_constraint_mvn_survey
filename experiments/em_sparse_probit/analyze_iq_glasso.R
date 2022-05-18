library(here)
library(data.table)

# 3d problem ------------------------------------------------------------

result_path = here("experiments", "em_sparse_probit", "output", "3d")

# preprocess --

method_paths = dir(result_path, full.names = TRUE)

prec_est_results = vector(mode="list", length = length(method_paths))
glasso_llik_results = vector(mode="list", length = length(method_paths))
                 
for (i in 1:length(method_paths)) {
  
  result_dir = dir(method_paths[i], full.names = TRUE)
  method_results = dir(method_paths[i], full.names = TRUE)
  
  prec_est_temp = vector(mode="list", length = length(method_results))
  glasso_llik_temp = vector(mode="list", length = length(method_results))
  
  for (j in 1:length(method_results)) {
    
    mr = method_results[[j]]
    fit = readRDS(mr)
    
    prec_est_temp[[j]] = data.table(
      "est_12" = fit$Prec[1, 2],
      "runtime" = attr(fit, "runtime"),
      "method" = attr(fit, "method"), 
      "rep" = j
    )
    
    glasso_llik_temp[[j]] = data.table(
      "iter" = 1:length(fit$gllaso_llik),
      "llik" = fit$gllaso_llik,
      "est_12" = sapply(fit$Prec_history, function(Prec) Prec[1, 2]),
      "method" = attr(fit, "method"),
      "rep" = j
    )
  }
  
  prec_est_results[[i]] = rbindlist(prec_est_temp)
  glasso_llik_results[[i]] = rbindlist(glasso_llik_temp)
}


# convergence check

glasso_llik_dt = rbindlist(glasso_llik_results)

library(ggplot2)

ggplot(glasso_llik_dt, aes(x = iter, y = llik, color = as.factor(rep))) +
  geom_line() +
  facet_wrap(vars(method)) +
  guides(color = "none")

# parameter estimate

prec_est_dt = rbindlist(prec_est_results)

ggplot(prec_est_dt, aes(x = method, y = est_12)) +
  geom_point()

# parameter estimate by iter
ggplot(glasso_llik_dt, aes(x = iter, y = est_12, color = as.factor(rep))) +
  geom_line() +
  facet_wrap(vars(method))

# why does llik look like the estimate?
ggplot(glasso_llik_dt, aes(x = llik, y = est_12, color = as.factor(rep))) +
  geom_point(pch =15) +
  facet_wrap(vars(method))
