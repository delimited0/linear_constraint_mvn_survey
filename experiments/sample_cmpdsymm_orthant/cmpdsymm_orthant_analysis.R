library(here)
library(data.table)
library(coda)
library(mcmcse)
library(jsonlite)

# test samples ------------------------------------------------------------
dims = read_json(
  here("experiments", "sample_cmpdsymm_orthant", "test_dim_conf.json"), 
  simplifyVector = TRUE
)

sample_path = here("experiments", "sample_cmpdsymm_orthant", "test_samples")

all_stats = rbindlist(lapply(
  dir(sample_path, full.names = TRUE), 
  function(method_path) {
    
    rbindlist(lapply(
      dir(method_path, full.names = TRUE), 
      function(dim_samples) {
        
        samples = readRDS(dim_samples)
        d = ncol(samples)
        
        point_est_summary = mcse.mat(samples)
        
        data.table(
          problem_d = d,
          d = 1:d,
          ess = effectiveSize(samples),
          mc_est = point_est_summary[, "est"],
          mc_se = point_est_summary[, "se"],
          runtime = attr(samples, "runtime"),
          method = attr(samples, "method")
        )
      }
    ))     
  }
))

# mean relative error (EP mean as reference)
ep_mean = rbindlist(lapply(dims, function(d) {
  mu = rep(0, d)
  Sigma = .5 * diag(d) + .5 * rep(1, d) %*% t(rep(1, d))
  lb = rep(0, d)
  ub = rep(Inf, d)
  
  moments = epmgpr::moments(lb, ub, mu, Sigma)
  
  data.table(
    mu = moments$mu[1],
    problem_d = d
  )
}))

all_stats = merge(all_stats, ep_mean)

avg_stats = all_stats[,
                      .(
                        mean_rel_error = mean(abs((mc_est - mu) / mu)),
                        ess_per_sec = mean(ess / runtime),
                        runtime = first(runtime)
                      ), 
                      by = list(method, problem_d)]

# visualization -----------------------------------------------------------
library(ggplot2)
library(ggrepel)

n_methods <- length(unique(avg_stats$method))
method_colors <- setNames(viridis::viridis_pal(option = "D")(n_methods), 
                          nm = unique(avg_stats$method))

perf_dim_style = list(
  theme_bw(),
  geom_text_repel(aes(label = label), color = "black"),
  guides(fill="none", color="none"),
  scale_color_manual(values = method_colors)
)

avg_stats[, label := ifelse(problem_d == max(problem_d), method, NA_character_),
          by = method]

# mean estimation accuracy
ggplot(avg_stats, aes(x = problem_d, y = mean_rel_error, 
                      fill = method, color = method)) +
  geom_point() + geom_line() +
  perf_dim_style

# mean effective samples / second
ggplot(avg_stats, aes(x = problem_d, y = log(ess_per_sec),
                      fill = method, color = method)) +
  geom_point() + geom_line() +
  perf_dim_style

# runtime
ggplot(avg_stats, aes(x = problem_d, y = runtime,
                      fill = method, color = method)) +
  geom_point() + geom_line() +
  perf_dim_style





