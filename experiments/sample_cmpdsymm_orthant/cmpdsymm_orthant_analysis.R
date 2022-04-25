library(here)
library(data.table)
library(coda)
library(mcmcse)
library(jsonlite)

# pick samples ------------------------------------------------------------
# dims = read_json(
#   here("experiments", "sample_cmpdsymm_orthant", "test_dim_conf.json"), 
#   simplifyVector = TRUE
# )
# sample_path = here("experiments", "sample_cmpdsymm_orthant", "test_samples")

dims = read_json(
  here("experiments", "sample_cmpdsymm_orthant", "dim_conf.json"), 
  simplifyVector = TRUE
)
sample_path = here("experiments", "sample_cmpdsymm_orthant", "samples")


# read samples --------------------------------------------------------------

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

# preprocess --------------------------------------------------------------

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

# save stats --------------------------------------------------------------
all_stat_path = 
  here("experiments", "sample_cmpdsymm_orthant", "samples", "all_stats.RDS")
saveRDS(all_stats, all_stat_path)

all_stats = readRDS(all_stat_path)

# visualization -----------------------------------------------------------
library(ggplot2)
library(ggrepel)
library(shades)

avg_stats = all_stats[,
                      .(
                        mean_rel_error = mean(abs((mc_est - mu) / mu)),
                        ess_per_sec = mean(ess / runtime),
                        runtime = first(runtime)
                      ), 
                      by = list(method, problem_d)]

n_methods <- length(unique(avg_stats$method))

method_colors <- setNames(viridis::viridis_pal(option = "D")(n_methods), 
                          nm = unique(avg_stats$method))
method_rgb = col2rgb(method_colors, )
label_colors = rgb(t(method_rgb / 2), maxColorValue = 255, 
                   names = colnames(method_rgb))

method_shapes = (0:14)[1:n_methods]
names(method_shapes) = names(method_colors)

perf_dim_style = list(
  theme_bw(),
  geom_text_repel(aes(label = label), force = 2),
  guides(fill="none", color="none", shape="none"),
  scale_color_manual(values = brightness(method_colors, delta(-.1))),
  scale_shape_manual(values = method_shapes)
)

avg_stats[, label := ifelse(problem_d == max(problem_d), method, NA_character_),
          by = method]

# mean estimation accuracy
ggplot(avg_stats, aes(x = problem_d, y = mean_rel_error, 
                      fill = method, color = method, shape = method)) +
  geom_point(size=2) + geom_line(linetype=2) +
  perf_dim_style +
  scale_y_log10() +
  labs(x = "Problem dimension", y = "Mean marginal relative error")

# mean effective samples / second
ggplot(avg_stats, aes(x = problem_d, y = ess_per_sec,
                      fill = method, color = method, shape = method)) +
  geom_point() + geom_line() +
  perf_dim_style +
  scale_y_log10() +
  labs(x = "Problem dimension", y = "Effective samples per second")

# runtime
ggplot(avg_stats, aes(x = problem_d, y = runtime,
                      fill = method, color = method)) +
  geom_point() + geom_line() +
  perf_dim_style +
  labs(x = "Problem dimension", y = "Seconds")





