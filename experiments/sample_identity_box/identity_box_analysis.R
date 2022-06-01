library(here)
library(data.table)
library(coda)
library(mcmcse)

marginal_ks <- function(samples, mu, Sigma, lb, ub) {
  d <- length(mu)
  ks <- rep(NA, d)
  for (k in 1:d) {
    m <- mu[k]
    s <- sqrt(Sigma[k, k])
    a <- lb[k]
    b <- ub[k]
    alp = (a - m) / s
    bet = (b - m) / s
    marginal_mean <-
      mu + m + s * (dnorm(alp) - dnorm(bet)) / (pnorm(bet) - pnorm(alp)) 
    marginal_sd <- 
      s^2 * (
        1 + ((alp * dnorm(alp) - bet * dnorm(bet)) / (pnorm(bet) - pnorm(alp))) -
          ((dnorm(alp) - dnorm(bet)) / (pnorm(bet) - pnorm(alp)))^2
      )
    ks[k] <- ks.test(samples[, k], truncnorm::ptruncnorm, 
                     a = lb[k], b = ub[k], mean = marginal_mean, 
                     sd = marginal_sd)$statistic
    ks[k] <- ks.test(samples[, k], truncnorm::ptruncnorm,
                     a = lb[k], b = ub[k], 
                     mean = mu[k], sd = sqrt(Sigma[k, k]))$statistic
  }
  return(ks)
}


# test sampling -----------------------------------------------------------
# sample_path = 
#   here("experiments", "sample_identity_box", "test_samples", "width=3")
sample_path = 
  here("experiments", "sample_identity_box", "samples", "width=3")


# preprocess --------------------------------------------------------------

# very slow... effectiveSize's fault I think
all_stats = rbindlist(lapply( 
  dir(sample_path, full.names = TRUE), 
  function(method) {
    
    method_stats = rbindlist(lapply(
      dir(method, full.names = TRUE), 
      function(dim_samples) {
        
        print(dim_samples)
        
        samples = readRDS(dim_samples)
        d = ncol(samples)
        half_width = attr(samples, "half_width")
        
        data.table(
          ks = marginal_ks(samples, 
                           mu = rep(0, d), Sigma = diag(d), 
                           lb = rep(-half_width, d), ub = rep(half_width, d)),
          ess = effectiveSize(samples),
          runtime = attr(samples, "runtime"),
          method = attr(samples, "method"),
          half_width = attr(samples, "half_width"),
          problem_d = d,
          d = 1:d
        )
      }
    ))
    
    method_stats
  }
))


# save stats --------------------------------------------------------------
all_stat_path = 
  here("experiments", "sample_identity_box", "samples", "all_stats.RDS")
saveRDS(all_stats, all_stat_path)

all_stats = readRDS(all_stat_path)

# visualization -----------------------------------------------------------
library(ggplot2)
library(ggrepel)
library(shades)


families = fread(here("sampling_method_directory.csv"))
all_stats = merge(all_stats, families)

n_methods <- length(unique(all_stats$method))
# method_colors <- setNames(viridis::viridis_pal(option = "D")(n_methods), 
#                           nm = unique(all_stats$method))
method_colors = setNames(RColorBrewer::brewer.pal(n_methods, "BrBG"), 
                         nm = unique(all_stats$method))

# methods by shape
shape_set = c(15:18, 7:14)
method_shapes = setNames(shape_set[1:n_methods], nm = unique(all_stats$method))

# family by color
n_family = length(unique(all_stats$family))
family_colors = setNames(viridis::turbo(n_family),
                         nm = unique(all_stats$family))

# common plot style
perf_dim_style = list(
  theme_bw(base_size = 14),
  geom_text_repel(aes(label = label), 
                  # color = "black", 
                  max.overlaps = Inf,
                  force = 10, 
                  # xlim = c(0, 3000),
                  nudge_x = 200),
  guides(fill="none", 
         # color="none",
         shape="none", linetype="none", color = "none"),
  xlim(0, 2500),
  # scale_color_manual(values = method_colors)
  scale_color_manual(values = family_colors),
  scale_shape_manual(values = method_shapes),
  theme(legend.position = "bottom")
)

avg_marginal_stats = all_stats[, 
                               .(ks = mean(ks), 
                                 ess = mean(ess), 
                                 runtime = min(runtime)),
                               by = list(method, family, problem_d)]

avg_marginal_stats[, 
                   label := ifelse(problem_d == max(problem_d), 
                                   method, NA_character_),
                   by = method]

# runtime ----
runtime_plot = ggplot(avg_marginal_stats, aes(x = problem_d, y = runtime, 
                               shape = method,
                               color = family)) +
  geom_point(size=3) + geom_line(linetype=2) +
  scale_y_log10(labels = function(x) format(x, scientific=FALSE)) +
  perf_dim_style +
  labs(x = "Dimension", y = "Runtime (seconds)")

ggsave(
  runtime_plot,
  filename = "runtime_plot.pdf",
  device = "pdf",
  path =  here("plots", "sample_identity_box"),
  width = 5,
  height = 6
)

# average marginal effective samples per second by problem dimension ----
esspersec_plot = ggplot(avg_marginal_stats, 
       aes(x = problem_d, y = ess / runtime, 
           shape = method, color = family)) +
  scale_y_log10(labels = function(x) format(x, scientific=FALSE)) +
  geom_point(size=3) + geom_line(linetype=2) +
  perf_dim_style +
  labs(x = "Dimension", y = "Mean marginal effective samples per second")
  
ggsave(
  esspersec_plot,
  filename = "essps_plot.pdf",
  device = "pdf",
  path =  here("plots", "sample_identity_box"),
  width = 5,
  height = 7
)
  
# average marginal ks distance by problem dimension ----
ks_plot = ggplot(avg_marginal_stats, aes(x = problem_d, ks, 
                               shape = method, color = family)) +
  geom_point(size=3) + geom_line(linetype=2) +
  perf_dim_style +
  labs(x = "Dimension", y = "Mean marginal KS test statistic")

ggsave(
  ks_plot,
  filename = "ks_plot.pdf",
  device = "pdf",
  path =  here("plots", "sample_identity_box"),
  width = 5,
  height = 7
)

