####################################################
### Spatial model diagnostics
####################################################
### Meuse: Compute SPEP of RF model
####################################################
### (c) 2021 Alexander Brenning
####################################################

library("sperrorest")
library("purrr")
library("magrittr")
library("alexmisc")
library("future")
library("future.callr")

plan(future.callr::callr, workers = 10)

# Additional functions that are not included in packages:
source("spdiagnostics-functions.R", encoding = "UTF-8")

# Load data and model settings:
load("run_meuse.rda")

# Arguments for partition_discs2:
smp_args <- list(maxdist = MAXDIST, 
                 repetition = NREP, 
                 seed1 = 123)

res_lm <- sperrorest(fo, d, 
                     model_fun = lm, 
                     coords = c("utmx", "utmy"),
                     smp_fun = partition_discs2,
                     smp_args = smp_args,
                     importance = TRUE, imp_permutations = NPERM,
                     imp_sample_from = "all",
                     distance = TRUE, 
                     verbose = 1, progress = FALSE,
                     err_fun = err_meuse)

save(res_lm, file = gsub("_res", "_res_lm", OUTFILE))

#print(summary(res_lm$error_rep)["test_rmse", "mean"])
#print(summary(res_lm$importance)[, c("mean.rmse", "sd.rmse")])

cat("\nDone.\n")
