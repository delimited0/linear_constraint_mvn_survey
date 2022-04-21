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

n_methods <- length(unique(all_stats$method))
method_colors <- setNames(viridis::viridis_pal(option = "D")(n_methods), 
                          nm = unique(all_stats$method))

perf_dim_style = list(
  theme_bw(),
  geom_text_repel(aes(label = label), color = "black", max.overlaps = Inf),
  guides(fill="none", color="none"),
  scale_color_manual(values = method_colors)
)

avg_marginal_stats = all_stats[, 
                               .(ks = mean(ks), 
                                 ess = mean(ess), 
                                 runtime = min(runtime)),
                               by = list(method, problem_d)]

avg_marginal_stats[, 
                   label := ifelse(problem_d == max(problem_d), 
                                   method, NA_character_),
                   by = method]

# runtime ----
ggplot(avg_marginal_stats, aes(x = problem_d, y = runtime, color = method)) +
  geom_point() + geom_line() +
  perf_dim_style

# average marginal effective samples per second by problem dimension ----
ggplot(avg_marginal_stats, aes(x = problem_d, y = log(ess / runtime), color = method)) +
  geom_point() + geom_line() +
  perf_dim_style
  
# average marginal ks distance by problem dimension ----
ggplot(avg_marginal_stats, aes(x = problem_d, ks, color = method)) +
  geom_point() + geom_line() +
  perf_dim_style


