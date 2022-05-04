# must install renv first ...
install.packages("renv")

# figure out which packages are from CRAN, which from github
custom_pkgs = c(
  "sampletmvn",
  "lincongauss",
  "met",
  "epmgpr",
  "tmg",
  "fftwtools",
  "hccmvn"
)

deps = renv::dependencies()
cran_pkgs = setdiff(unique(deps$Package), custom_pkgs)

# install packages from CRAN necessary to run experiments
for (pkg in cran_pkgs) {
  
  find_path = find.package(pkg, quiet=TRUE)  
  
  if (length(find_path) == 0) {  # no such package
    install.packages(pkg)
  }
}

# renv::install(cran_pkgs)



