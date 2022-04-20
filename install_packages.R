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
renv::install(cran_pkgs)

# install packages no longer available on CRAN
renv::install("cran/tmg@bd996adcc584886cc66fc9cedd5d1426e54f63b1")

# install packages from Github necessary to run experiments
renv::install("delimited0/epmgpr")
renv::install("delimited0/sampletmvn")
renv::install("delimited0/met")
renv::install("delimited0/rcpp-lin-con-gauss")
renv::install("JCatwood/HCCMVN")


