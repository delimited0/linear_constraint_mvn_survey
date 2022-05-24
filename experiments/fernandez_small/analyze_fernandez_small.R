library(here)
library(data.table)

# pick results ------------------------------------------------------------

result_path = here("experiments", "fernandez_small", "results")
# result_path = here("experiments", "fernandez_small", "test_results")

# preprocess --------------------------------------------------------------

method_paths = dir(result_path, full.names = TRUE)
results = vector(mode="list", length = length(method_paths))

for (i in 1:length(method_paths)) {
  
  method_results = dir(method_paths[i], full.names = TRUE)
  
  results[[i]] = rbindlist(lapply(method_results, function(mr) {
    
    est = readRDS(mr)
    attributes(est) = NULL
    
    error = sapply(est, function(e) {
      err = attr(e, "error")
      if (is.null(err)) err = NA
      return(err)
    })
    
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

# relative error ----------------------------------------------------------
bounds =
  rbind(data.table(type = "Lower bound",
                   d = c(2, 3, 5, 10, 15, 20, 25, 30, 40, 50),
                   value = c(.0148955, .0010771, 2.4505 * 10^-6, 8.5483 * 10^-15,
                             1.3717 * 10^-25, 1.7736 * 10^-38, 2.674 * 10^-53, 
                             6.09*10^-70, 2.17 * 10^-108, 2.131 * 10^-153)),
        data.table(type = "Upper bound",
                   d = c(2, 3, 5, 10, 15, 20, 25, 30, 40, 50),
                   value = c(.0149, .00108, 2.48 * 10^-6, 2.1046 * 10^-14,
                             1.43 * 10^-25, 1.869 * 10^-38, 2.83 * 10^-53,
                             6.46 * 10^-70, 2.3 * 10^-108, 2.24 * 10^-153)))
bounds = dcast(bounds, d ~ type)

all_stats = merge(all_stats, bounds, by = "d")

families = fread(here("prob_method_directory.csv"))
all_stats = merge(all_stats, families, by = "method")

# how far out of the bounds are each methods' estimates?
all_stats[, c("ldist", "udist") := list(max(`Lower bound` - estimate, 0),
                                        max(estimate - `Upper bound`, 0)),
          by = c("method", "d", "rep")][
            , max_dist := max(ldist, udist), by = c("method", "d", "rep")]
all_stats[, relerror := max_dist / `Lower bound`]
all_stats[, label := ifelse(d == max(d), method, NA_character_), by = method]

avg_stats = all_stats[,
                      .(
                        estimate = mean(estimate),
                        runtime = mean(runtime),
                        relerror = mean(relerror),
                        se_relerror = sd(relerror) / sqrt(.N),
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
                  nudge_x = 2),
  guides(fill="none", 
         # color="none",
         shape="none", linetype="none"),
  scale_color_manual(values = family_colors),
  scale_shape_manual(values = method_shapes),
  theme(legend.position = "bottom")
)

# runtime ----
runtime_plot =
  ggplot(avg_stats, 
         aes(x = d, y = runtime, shape = method, color = family)) +
  geom_point(size = 2) + geom_line(linetype = 2) +
  geom_linerange(aes(ymin = runtime - 2*se_runtime, ymax = runtime + 2*se_runtime)) +
  perf_dim_style +
  scale_y_log10(labels = function(x) format(x, scientific=FALSE)) +
  xlim(0, 60) + 
  labs(x = "Dimension", y = "Runtime (seconds)")

ggsave(runtime_plot,
       filename = "runtime.pdf",
       device = "pdf",
       path = here("plots", "prob_fernandez_small"),
       width = 6,
       height = 6)

# accuracy ----
accuracy_plot = ggplot(avg_stats[!(method %in% c("uvcdn", "bvcdn", "dvcdn"))],
       aes(x = d, y = relerror, color = method)) +
  geom_point() + geom_line() +
  perf_dim_style +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Dimension", y = "Max relative bound exceedance")
  # scale_y_log10() +

ggsave(accuracy_plot,
       filename = "accuracy.pdf",
       device = "pdf",
       path = here("plots", "prob_fernandez_small"),
       width = 6,
       height = 6)

