library(here)
library(data.table)

# pick results ------------------------------------------------------------

result_path = here("experiments", "prob_exp_covariance", "results", "dep=0.1")

# preprocess --------------------------------------------------------------

method_paths = dir(result_path, full.names = TRUE)
results = vector(mode="list", length = length(method_paths))

for (i in 1:length(method_paths)) {
  
  method_results = dir(method_paths[i], full.names = TRUE)
  
  results[[i]] = rbindlist(lapply(method_results, function(mr) {
    
    # grab dim size from file name
    dim_spec = stringr::str_extract(mr, "d=\\d{1,}")
    d = as.numeric(
      stringr::str_extract(dim_spec, "\\d{1,}")
    )
    
    est = readRDS(mr)
    error = attr(est, "error")
    if (is.null(error))
      error = NA
    
    dt = data.table(
      "estimate" = est,
      "error" = error,
      "method" = attr(est, "method"),
      "runtime" = attr(est, "runtime"),
      # "d" = attr(est, "d")
      d = d
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

all_stats[, label := ifelse(d == max(d), method, NA_character_), by = method]

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
  xlim(0, 18000),
  scale_color_manual(values = family_colors),
  scale_shape_manual(values = method_shapes),
  theme(legend.position = "bottom")
)


# runtime -----
ggplot(all_stats, aes(x = d, y = runtime, shape = method, color = family)) +
  geom_point(size = 2) + geom_line(linetype = 2) +
  perf_dim_style +
  scale_y_log10()

# accuracy ----
ggplot(all_stats, aes(x = d, y = estimate, shape = method, color = family)) +
  geom_point(size = 2) + geom_line(linetype = 2) +
  geom_errorbar(aes(ymin = estimate - 2*error, ymax = estimate + 2*error)) +
  perf_dim_style
  




  
