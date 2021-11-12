####################################################
### Spatial model diagnostics
####################################################
### Meuse: Compute SPEP of OK-RF model
####################################################
### (c) 2021 Alexander Brenning
####################################################

library("sperrorest")
library("randomForest")
library("purrr")
library("magrittr")
library("alexmisc")
library("gstat")
library("future")
library("future.callr")

#plan(future.callr::callr, workers = 10)

# Additional functions that are not included in packages:
source("spdiagnostics-functions.R", encoding = "UTF-8")

# Load data and model settings:
load("run_meuse.rda")

# Arguments for partition_discs2:
smp_args <- list(maxdist = MAXDIST, 
                 repetition = NREP, 
                 seed1 = 123)

# fit <- okrf(fo, d, ntree = NTREE,
#             locations = ~utmx+utmy, 
#             range = 500,
#             taper_radius = 500)

res_orf <- sperrorest(fo, d, model_fun = okrf, 
                      model_args = list(ntree = NTREE,
                                        locations = ~utmx+utmy, 
                                        range = 500,
                                        taper_radius = 500),
                      coords = c("utmx", "utmy"),
                      smp_fun = partition_discs2,
                      smp_args = smp_args,
                      importance = TRUE, imp_permutations = NPERM,
                      imp_sample_from = "all",
                      distance = TRUE, 
                      verbose = 1, progress = FALSE,
                      mode_rep = "loop", mode_fold = "loop",
                      err_fun = err_meuse)

save(res_orf, file = gsub("_res", "_res_okrf", OUTFILE))

options(digits = 3)
print(summary(res_orf$error_rep))
print(summary(res_orf$importance))

cat("\nDone.\n")
