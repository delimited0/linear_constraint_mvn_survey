library(here)
library(data.table)

# pick results ------------------------------------------------------------

result_path = here("experiments", "fernandez_small", "test_results")

# preprocess --------------------------------------------------------------

method_paths = dir(result_path, full.names = TRUE)
results = vector(mode="list", length = length(method_paths))

for (i in 1:length(method_paths)) {
  
  method_results = dir(method_paths[i], full.names = TRUE)
  
  results[[i]] = rbindlist(lapply(method_results, function(mr) {
    
    est = readRDS(mr)
    error = attr(est, "error")
    if (is.null(error))
      error = NA
    
    dt = data.table(
      "estimate" = est,
      "error" = error,
      "method" = attr(est, "method"),
      "runtime" = attr(est, "runtime"),
      "d" = attr(est, "d")
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

# how far out of the bounds are each methods' estimates?
all_stats[, c("ldist", "udist") := list(max(`Lower bound` - estimate, 0),
                                        max(estimate - `Upper bound`, 0)),
          by = c("method", "d")][
            , max_dist := max(ldist, udist), by = c("method", "d")]
all_stats[, relerror := max_dist / `Lower bound`]
all_stats[, label := ifelse(d == max(d), method, NA_character_), by = method]

ggplot(all_stats[!(method %in% c("uvcdn", "bvcdn", "dvcdn"))],
       aes(x = d, y = relerror, color = method)) +
  geom_point() + geom_line() +
  geom_text_repel(aes(label = label), na.rm = TRUE, show.legend = FALSE) +
  # scale_y_log10() +
  theme_bw() + 
  guides(fill=FALSE, color=FALSE)
