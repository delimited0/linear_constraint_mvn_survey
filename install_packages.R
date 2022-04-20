# must install renv first ...
install.packages("renv")

# figure out which packages are from CRAN, which from github
custom_pkgs = c(
  "sampletmvn",
  "lincongauss",
  "met",
  "epmgpr"
)

deps = renv::dependencies()
cran_pkgs = setdiff(unique(deps$Package), custom_pkgs)

# install packages from CRAN necessary to run experiments
install.packages(cran_pkgs)

# install packages from Github necessary to run experiments
devtools::install_github("")