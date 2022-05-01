library(here)
library(data.table)

# pick results ------------------------------------------------------------

# correlation = 0.5
result_path = here("experiments", "prob_cmpdsymm_orthant", "results", "corr=0.5")


# preprocess --------------------------------------------------------------

method_paths = dir(result_path, full.names = TRUE)

results = vector(mode="list", length = length(method_paths))

for (i in 1:length(method_paths)) {
  
  result_dir = dir(method_paths[i], full.names = TRUE)
  
  results[[i]] = rbindlist(lapply(
    result_dir, function(res) {
      df = readRDS(result_dir)
      df
    }
  ))
}

sapply(corr_dirs, function(path) dir(path, full.names=TRUE))

lapply(
  dir(result_path, full.names=TRUE),
  function(method_path) {
    
    lapply(
      dir(method_path, full.names = TRUE)
  }

  
  dir(dir(dir(result_path, full.names=TRUE), full.names=TRUE)[1]kjj