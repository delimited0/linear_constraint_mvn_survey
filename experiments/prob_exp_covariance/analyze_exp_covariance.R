library(here)
library(data.table)

# pick results ------------------------------------------------------------

# result_path = here("experiments", "prob_exp_covariance", "results", "dep=0.1")
result_path = here("experiments", "prob_exp_covariance", "test_results", "dep=0.1")

# preprocess --------------------------------------------------------------

method_paths = dir(result_path, full.names = TRUE)
results = vector(mode="list", length = length(method_paths))

for (i in 1:length(method_paths)) {
  
  method_results = dir(method_paths[i], full.names = TRUE)
  
  results[[i]] = rbindlist(lapply(method_results, function(mr) {
    
    # grab dim size from file name
    # dim_spec = stringr::str_extract(mr, "d=\\d{1,}")
    # d = as.numeric(
    #   stringr::str_extract(dim_spec, "\\d{1,}")
    # )
    # 
    est = readRDS(mr)
    attributes(est) = NULL
    
    error = sapply(est, function(e) {
      err = attr(e, "error")
      if (is.null(err)) err = NA
      return(err)
    })
    if (is.null(error))
      error = NA
    
    dt = data.table(
      "estimate" = c(est, recursive=TRUE),
      "error" = error,
      "method" = sapply(est, function(e) attr(e, "method")),
      "runtime" = sapply(est, function(e) attr(e, "runtime")),
      "d" = sapply(est, function(e) attr(e, "d")),
      "rep" = sapply(est, function(e) attr(e, "rep"))
      # d = d
    )  
    
    return(dt)
  }), fill = TRUE)
}

all_stats = rbindlist(results)


# visualization -----------------------------------------------------------
library(ggplot2)
library(ggrepel)

families = fread(here("prob_method_directory.csv"))
all_stats = merge(all_stats, families, by = "method")

avg_stats = all_stats[,
                      .(
                        estimate = mean(estimate),
                        runtime = mean(runtime),
                        se_estimate = sd(estimate) / sqrt(.N),
                        se_runtime = sd(runtime) / sqrt(.N)
                      ),
                      by = list(method, d, family)]

avg_stats[, label := ifelse(d == max(d), method, NA_character_), by = method]

# methods by shape
n_methods = length(unique(all_stats$method))
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
  xlim(0, 1800),
  scale_color_manual(values = family_colors),
  scale_shape_manual(values = method_shapes),
  theme(legend.position = "bottom")
)

# runtime -----
runtime_plot = ggplot(avg_stats, 
                      aes(x = d, y = runtime, shape = method, color = family)) +
  geom_point(size = 2) + geom_line(linetype = 2) +
  geom_linerange(aes(ymin = runtime - 2*se_runtime, ymax = runtime + 2*se_runtime)) +
  perf_dim_style +
  scale_y_log10(labels = function(x) format(x, scientific=FALSE)) +
  labs(x = "Dimension", y = "Runtime (seconds)")

ggsave(runtime_plot,
       filename = "runtime.pdf",
       device = "pdf",
       path = here("plots", "prob_exp_covariance"),
       width = 6,
       height = 6)

# accuracy ----
accuracy_plot = ggplot(avg_stats, 
                       aes(x = d, y = estimate, 
                           shape = method, color = family)) +
  geom_point(size = 2) + geom_line(linetype = 2) +
  geom_linerange(aes(ymin = estimate - 2*se_estimate, 
                     ymax = estimate + 2*se_estimate)) +
  perf_dim_style +
  labs(x = "Dimension", y = "Estimate")

ggsave(accuracy_plot,
       filename = "accuracy.pdf",
       device = "pdf",
       path = here("plots", "prob_exp_covariance"),
       width = 6,
       height = 6)
  




  
