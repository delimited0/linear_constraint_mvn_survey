library(jsonlite)
library(here)

target_widths = function(dims, target_prob) {
  
  f = function(y) pnorm(y) - pnorm(-y) - target_prob^(1/d)
  
  widths = rep(NA, length(dims))
  for (i in 1:length(dims)) {
    d = dims[i]
    opt = uniroot(f, c(-10, 10))
    widths[i] = opt$root
  }
  
  return(widths)
} 


dims = c(10, 30, 50, 100, 300, 500, 1000, 2000, 3000, 4000)
tw = target_widths(dims, .5)

write_json(
  x = list(dims = dims, 
           half_widths = tw),
  path = here("experiments", "sample_identity_box", "foo.json")
)

dims = c(2, 5, 10, 50)
tw = target_widths(dims, .5)
