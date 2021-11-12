####################################################
### Spatial model diagnostics
####################################################
### Install required packages from CRAN and Github
####################################################
### (c) 2021 Alexander Brenning
####################################################

### Install required packages from CRAN:

pckgs <- c(
  "devtools",      # package management
  "sp",            # Meuse data
  "raster",        # raster handling
  "gstat",         # kriging
  "randomForest",  # random forest
  "spgwr",         # GWR
  "dplyr",         # data handling
  "purrr",         # data handling
  "magrittr",      # data handling
  "future",        # parallelization
  "future.callr"   # parallelization
)

for (pckg in pckgs) {
  ipckgs <- rownames(installed.packages())
  if (all(ipckgs != pckg)) {
    try(install.packages(pckg))
  }
}


### Install additional required packages / most recent
### development versions from Github:

library("devtools")

# Interpretable machine learning:
devtools::install_github("christophM/iml")
devtools::install_github("alexanderbrenning/wiml")

# Spatial cross-validation:
devtools::install_github("giscience-fsu/sperrorest")

# Miscellaneous functions:
devtools::install_github("alexanderbrenning/alexmisc")
