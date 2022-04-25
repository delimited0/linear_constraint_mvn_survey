library(here)
library(data.table)
library(coda)
library(mcmcse)

# pick samples ------------------------------------------------------------

# sample_path = here("experiments", "proton_radius", "test_samples")
sample_path = here("experiments", "proton_radius", "samples")


# preprocess --------------------------------------------------------------

all_stats = rbindlist(lapply(
  dir(sample_path, full.names = TRUE), 
  function(method_path) {
    
    rbindlist(lapply(
      dir(method_path, full.names = TRUE), 
      function(method_rep) {
        
        samples = readRDS(method_rep)
        d = ncol(samples)
        
        # how does Q_max enter here?
        rp_samples = sqrt(-6 * samples[-c(1:100), 2]) 
        
        point_est_summary = mcse(rp_samples)
        
        data.table(
          posterior_mean = point_est_summary$est,
          posterior_median = median(rp_samples),
          mc_se = point_est_summary$se,
          runtime = attr(samples, "runtime"),
          method = attr(samples, "method"),
          ess = effectiveSize(rp_samples)
        )
      }
    ))     
  }
))


# visualization -----------------------------------------------------------
library(ggplot2)

# true parameter recovery
recovery_plot = ggplot(all_stats[method != "epess"], 
                       aes(x = method, y = posterior_median)) +
  geom_point() + 
  geom_hline(yintercept = 0.84) +
  theme_bw() +
  guides(x = guide_axis(angle = 45)) +
  labs(x = "Method", y = "Posterior median")

ggsave(
  recovery_plot,
  filename = "recovery_plot.pdf",
  device = "pdf",
  path =  here("plots", "proton_radius"),
  width = 6,
  height = 4
)

# MC standard error
mcse_plot = ggplot(all_stats, aes(x = method, y = mc_se)) +
  geom_point() + 
  theme_bw() +
  guides(x = guide_axis(angle = 45)) +
  labs(x = "Method", y = "Monte Carlo standard error")

ggsave(
  mcse_plot,
  filename = "mcse_plot.pdf",
  device = "pdf",
  path =  here("plots", "proton_radius"),
  width = 6,
  height = 4
)

# ESS / second
esspersec_plot = ggplot(all_stats, aes(x = method, y = ess / runtime)) +
  geom_point() +
  theme_bw() +
  guides(x = guide_axis(angle = 45)) +
  labs(x = "Method", y = "Effective samples per second")

ggsave(
  esspersec_plot,
  filename = "essps_plot.pdf",
  device = "pdf",
  path =  here("plots", "proton_radius"),
  width = 6,
  height = 4
)

# runtime
runtime_plot = ggplot(all_stats, aes(x = method, y = runtime)) +
  geom_point() +
  theme_bw() +
  guides(x = guide_axis(angle = 90)) +
  labs(x = "Method", y = "Runtime (seconds)")

ggsave(
  runtime_plot,
  filename = "runtime_plot.pdf",
  device = "pdf",
  path =  here("plots", "proton_radius"),
  width = 6,
  height = 4
)




