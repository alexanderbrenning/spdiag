####################################################
### Spatial model diagnostics
####################################################
### Meuse: Compute SPEP of RF model
####################################################
### (c) 2021 Alexander Brenning
####################################################

library("sperrorest")
library("sp")
library("spgwr")
library("purrr")
library("magrittr")
library("alexmisc")
library("future")
library("future.callr")

# Additional functions that are not included in packages:
source("spdiagnostics-functions.R", encoding = "UTF-8")

# Load data and model settings:
load("run_meuse.rda")

# Arguments for partition_discs2:
smp_args <- list(maxdist = MAXDIST, 
                 repetition = NREP, 
                 seed1 = 123)

res_gwr <- sperrorest(fo, d, 
                      model_fun = mygwr, 
                      model_args = list(coords = ~utmx+utmy),
                      coords = c("utmx", "utmy"),
                      smp_fun = partition_discs2,
                      smp_args = smp_args,
                      importance = TRUE, imp_permutations = NPERM,
                      imp_sample_from = "all",
                      # do not run this in parallel mode:
                      mode_rep = "loop", mode_fold = "loop",
                      distance = TRUE, 
                      progress = FALSE, verbose = 0,
                      err_fun = err_meuse)

save(res_gwr, file = gsub("_res", "_res_gwr", OUTFILE))

cat("\nDone.\n")
