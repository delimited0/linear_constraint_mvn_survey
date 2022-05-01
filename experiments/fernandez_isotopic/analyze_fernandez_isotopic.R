library(here)
library(data.table)

# pick results ------------------------------------------------------------

result_path = here("experiments", "fernandez_isotopic", "test_results")


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
  }))
}

all_stats = rbindlist(results)


# visualization -----------------------------------------------------------
library(ggplot2)
library(ggrepel)


# relative error ----------------------------------------------------------
bounds <- 
  rbind(data.table(type = "Lower bound",
                   d = c(2, 3, 10, 20, 25, 50, 80, 100, 120, 150, 200, 250),
                   value = c(.09114, .02303, 1.338 * 10^-6, 1.080 * 10^-12,
                             9.77 * 10^-16, 5.925 * 10^-31, 3.252 * 10^-49,
                             2.18 * 10^-61, 1.462 * 10^-73, 8.026 * 10^-92,
                             2.954 * 10^-122, 1.087 * 10^-152)),
        data.table(type = "Upper bound", 
                   d = c(2, 3, 10, 20, 25, 50, 80, 100, 120, 150, 200, 250),
                   value = c(.09205, .0234, 1.454 * 10^-6, 1.289 * 10^-12,
                             1.222 * 10^-15, 9.368 * 10^-31, 6.812 * 10^-49,
                             5.5 * 10^-61, 4.45 * 10^-73, 3.23 * 10^-91,
                             1.905 * 10^-121, 1.12 * 10^-151)))
bounds <- dcast(bounds, d ~ type)

all_stats = merge(all_stats, bounds, by = "d")

# how far out of the bounds are each methods' estimates?
all_stats[, c("ldist", "udist") := list(max(`Lower bound` - estimate, 0),
                                        max(estimate - `Upper bound`, 0)),
          by = c("method", "d")][
            , max_dist := max(ldist, udist), by = c("method", "d")]
all_stats[, relerror := max_dist / `Lower bound`]
all_stats[, label := ifelse(d == max(d), method, NA_character_), by = method]

ggplot(all_stats, aes(x = d, y = relerror, color = method)) +
  geom_point() + geom_line() +
  scale_y_log10()
