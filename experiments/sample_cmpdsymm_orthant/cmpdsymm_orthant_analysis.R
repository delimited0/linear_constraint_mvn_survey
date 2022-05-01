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

families = fread(here("sampling_method_directory.csv"))
all_stats = merge(all_stats, families, by = "method")

avg_stats = all_stats[,
                      .(
                        mean_rel_error = mean(abs((mc_est - mu) / mu)),
                        ess_per_sec = mean(ess / runtime),
                        runtime = first(runtime)
                      ), 
                      by = list(method, problem_d, family)]

n_methods <- length(unique(avg_stats$method))

# methods by shape
shape_set = c(15:18, 7:14)
method_shapes = setNames(shape_set[1:n_methods], nm = unique(all_stats$method))

# family by color
n_family = length(unique(all_stats$family))
family_colors = setNames(viridis::turbo(n_family),
                         nm = unique(all_stats$family))

# common plot style
perf_dim_style = list(
  theme_bw(),
  geom_text_repel(aes(label = label), 
                  # color = "black", 
                  max.overlaps = Inf,
                  force = 10, 
                  # xlim = c(0, 3000),
                  nudge_x = 200),
  guides(fill="none", 
         # color="none",
         shape="none", linetype="none"),
  xlim(0, 4500),
  # scale_color_manual(values = method_colors)
  scale_color_manual(values = family_colors),
  scale_shape_manual(values = method_shapes),
  theme(legend.position = "bottom")
)

avg_stats[, label := ifelse(problem_d == max(problem_d), method, NA_character_),
          by = method]

# mean estimation accuracy ----
mean_est_plot = ggplot(avg_stats, aes(x = problem_d, y = mean_rel_error, 
                      color = family, shape = method)) +
  geom_point(size=3) + geom_line(linetype=2) +
  perf_dim_style +
  scale_y_log10() +
  labs(x = "Problem dimension", y = "Mean marginal relative error")

ggsave(
  mean_est_plot,
  filename = "mean_est_plot.pdf",
  device = "pdf",
  path =  here("plots", "sample_cmpdsymm_orthant"),
  width = 6,
  height = 5
)

# mean effective samples / second
esspersec_plot = ggplot(avg_stats, aes(x = problem_d, y = ess_per_sec,
                      color = family, shape = method)) +
  geom_point(size=3) + geom_line(linetype=2) +
  perf_dim_style +
  scale_y_log10() +
  labs(x = "Problem dimension", y = "Effective samples per second")

ggsave(
  esspersec_plot,
  filename = "essps_plot.pdf",
  device = "pdf",
  path =  here("plots", "sample_cmpdsymm_orthant"),
  width = 6,
  height = 5
)

# runtime
runtime_plot = ggplot(avg_stats, aes(x = problem_d, y = runtime,
                      color = family, shape = method)) +
  geom_point(size=3) + geom_line(linetype=2) +
  perf_dim_style +
  scale_y_log10() +
  labs(x = "Problem dimension", y = "Seconds")

ggsave(
  runtime_plot,
  filename = "runtime_plot.pdf",
  device = "pdf",
  path =  here("plots", "sample_cmpdsymm_orthant"),
  width = 6,
  height = 5
)





