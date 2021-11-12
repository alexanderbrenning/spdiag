####################################################
### Spatial model diagnostics
####################################################
### Meuse: Compute SPEP of OK model
####################################################
### (c) 2021 Alexander Brenning
####################################################

library("sperrorest")
library("gstat")
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

res_ok <- sperrorest(fo, d, 
                     model_fun = mykrigeOK,
                     model_args = list(range = 500, 
                                       locations = ~ utmx + utmy),
                     pred_args = list(nmax = NMAX),
                     coords = c("utmx", "utmy"),
                     smp_fun = partition_discs2,
                     smp_args = smp_args,
                     importance = TRUE, imp_permutations = NPERM,
                     imp_sample_from = "all",
                     distance = TRUE, 
                     progress = FALSE, verbose = 1, 
                     # do not run this in parallel mode:
                     mode_rep = "loop", mode_fold = "loop",
                     err_fun = err_meuse)

save(res_ok, file = gsub("_res", "_res_ok", OUTFILE))

cat("\nDone.\n")
