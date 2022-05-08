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

sapply(corr_dirs, function(path) dir(path, full.names=TRUE))

lapply(
  dir(result_path, full.names=TRUE),
  function(method_path) {
    
    lapply(
      dir(method_path, full.names = TRUE)
  }

  
  dir(dir(dir(result_path, full.names=TRUE), full.names=TRUE)[1]kjj