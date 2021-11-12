####################################################
### Spatial model diagnostics
####################################################
### Meuse: Compute SPEP of RF model
####################################################
### (c) 2021 Alexander Brenning
####################################################

library("sperrorest")
library("randomForest")
library("purrr")
library("magrittr")
library("alexmisc")
library("future")
library("future.callr")

plan(future.callr::callr, workers = 10)
# Workstation will run out of memory when using 20 workers...

# Additional functions that are not included in packages:
source("spdiagnostics-functions.R", encoding = "UTF-8")

# Load data and model settings:
load("run_meuse.rda")

# Arguments for partition_discs2:
smp_args <- list(maxdist = MAXDIST, 
                 repetition = NREP, 
                 seed1 = 123)

res_rf <- sperrorest(fo, d, 
                     model_fun = randomForest, 
                     model_args = list(ntree = NTREE),
                     coords = c("utmx", "utmy"),
                     smp_fun = partition_discs2,
                     smp_args = smp_args,
                     importance = TRUE, imp_permutations = NPERM,
                     imp_sample_from = "all",
                     distance = TRUE,
                     # mode_rep = "future", mode_fold = "sequential",
                     progress = FALSE, verbose = 1,
                     err_fun = err_meuse)

save(res_rf, file = gsub("_res", "_res_rf", OUTFILE))

# options(digits = 3)
# print(summary(res_rf$error_rep)["test_bias", "mean"])
# print(summary(res_rf$importance))

cat("\nDone.\n")
