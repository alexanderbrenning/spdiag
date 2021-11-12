####################################################
### Spatial model diagnostics
####################################################
### Meuse: Compute SPEP of NN model
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

# Additional functions that are not included in packages:
source("spdiagnostics-functions.R", encoding = "UTF-8")

# Load data and model settings:
load("run_meuse.rda")

# Arguments for partition_discs2:
smp_args <- list(maxdist = MAXDIST, 
                 repetition = NREP, 
                 seed1 = 123)

res_nn <- sperrorest(fo, d, model_fun = nn, 
                     coords = c("utmx", "utmy"),
                     smp_fun = partition_discs2,
                     smp_args = smp_args,
                     importance = TRUE, imp_permutations = NPERM,
                     imp_sample_from = "all",
                     distance = TRUE, 
                     # do not run this in parallel mode:
                     mode_rep = "loop", mode_fold = "loop",
                     progress = FALSE, verbose = 1,
                     err_fun = err_meuse)

save(res_nn, file = gsub("_res", "_res_nn", OUTFILE))

options(digits = 3)
print(summary(res_nn$error_rep))
#print(summary(res_nn$error_fold))
print(summary(res_nn$importance))

cat("\nDone.\n")
