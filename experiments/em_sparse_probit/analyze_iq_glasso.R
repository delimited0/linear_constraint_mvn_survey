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
      "(1,2)" = fit$Prec[1, 2],
      "(1,3)" = fit$Prec[1, 3],
      "(2,3)" = fit$Prec[2, 3],
      "runtime" = attr(fit, "runtime"),
      "method" = attr(fit, "method"), 
      "rep" = j
    )
    
    glasso_llik_temp[[j]] = data.table(
      "iter" = 1:length(fit$glasso_llik),
      "llik" = fit$glasso_llik,
      "(1,2)" = sapply(fit$Prec_history, function(Prec) Prec[1, 2]),
      "method" = attr(fit, "method"),
      "rep" = j
    )
  }
  
  prec_est_results[[i]] = rbindlist(prec_est_temp)
  glasso_llik_results[[i]] = rbindlist(glasso_llik_temp)
}


### visualization ###
library(ggplot2)

## convergence check
glasso_llik_dt = rbindlist(glasso_llik_results)

ggplot(glasso_llik_dt, aes(x = iter, y = llik, color = as.factor(rep))) +
  geom_line() +
  facet_wrap(vars(method)) +
  guides(color = "none")

## parameter estimate
prec_est_dt = rbindlist(prec_est_results)
prec_est_dt = melt(prec_est_dt, id.vars = c("method", "rep", "runtime"))
prec_est_dt[, true_value := ifelse(variable == "(1,2)", .5, 0)]

ggplot(prec_est_dt, aes(x = method, y = value)) +
  # geom_point() +
  geom_boxplot() +
  geom_hline(aes(yintercept = true_value)) + 
  facet_wrap(vars(variable)) +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(x = "Method", y = "Estimated precision")

## runtime 
ggplot(prec_est_dt, aes(x = method, y = runtime)) +
  # geom_point() +
  geom_boxplot() +
  facet_wrap(vars(variable)) +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = 45))
  
# parameter estimate by iter
ggplot(glasso_llik_dt, aes(x = iter, y = est_12, color = as.factor(rep))) +
  geom_line() +
  facet_wrap(vars(method))

# why does llik look like the estimate?
ggplot(glasso_llik_dt, aes(x = llik, y = est_12, color = as.factor(rep))) +
  geom_point(pch =15) +
  facet_wrap(vars(method))


# 5d sparse problem -------------------------------------------------------

result_path = here("experiments", "em_sparse_probit", "output", "5d_a")

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
    
    upper_tri_mask = upper.tri(fit$Prec)
    upper_tri = fit$Prec[upper_tri_mask]
    
    idx = which(upper_tri_mask, arr.ind=T) 
    idx_names = paste0(
      "(", 
      apply(idx, 1, function(x) paste0(x, collapse = ",")),
      ")"
    )
    
    est_dt = data.table(
      "runtime" = attr(fit, "runtime"),
      "method" = attr(fit, "method"), 
      "rep" = j
    )
    est_dt[, (idx_names) := as.list(upper_tri)]
        
    prec_est_temp[[j]] = est_dt
    
    glasso_llik_temp[[j]] = data.table(
      "iter" = 1:length(fit$glasso_llik),
      "llik" = fit$glasso_llik,
      "(3, 4)" = sapply(fit$Prec_history, function(Prec) Prec[3, 4]),
      "method" = attr(fit, "method"),
      "rep" = j
    )
  }
  
  prec_est_results[[i]] = rbindlist(prec_est_temp)
  glasso_llik_results[[i]] = rbindlist(glasso_llik_temp)
}


# convergence check
library(ggplot2)

glasso_llik_dt = rbindlist(glasso_llik_results)

ggplot(glasso_llik_dt, aes(x = iter, y = llik, color = as.factor(rep))) +
  geom_line() +
  facet_wrap(vars(method), scales = "free") +
  guides(color = "none")

# parameter estimate
prec_est_dt = rbindlist(prec_est_results)
prec_est_dt = melt(prec_est_dt, id.vars = c("method", "rep", "runtime"))

ggplot(prec_est_dt, aes(x = method, y = value)) +
  geom_point() +
  facet_wrap(vars(variable))
