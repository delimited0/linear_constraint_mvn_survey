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


# problem parameters ------------------------------------------------------


# test sampling -----------------------------------------------------------
sample_path = 
  here("experiments", "sample_identity_box", "test_samples", "width=3")

all_stats = rbindlist(lapply( 
  dir(sample_path, full.names = TRUE), 
  function(method) {
    
    method_stats = rbindlist(lapply(
      dir(method, full.names = TRUE), 
      function(dim_samples) {
        
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

# visualization -----------------------------------------------------------
library(ggplot2)

avg_marginal_stats = all_stats[, 
                               .(ks = mean(ks), ess = mean(ess)),
                               by = list(method, problem_d)]

# runtime
ggplot(all_stats, aes(x = problem_d, y = runtime, color = method)) +
  geom_point() + geom_line()

# average marginal ks distance by dimension
ggplot(avg_marginal_stats, aes(x = problem_d, ks, color = method)) +
  geom_point() + geom_line()


