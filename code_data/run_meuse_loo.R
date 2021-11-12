####################################################
### Spatial model diagnostics: Meuse
####################################################
### LOO-CV
####################################################
### (c) 2021 Alexander Brenning
####################################################

# Data handling:
library("purrr")
library("magrittr")

# Cross-validation:
library("sperrorest")

# Models:
library("randomForest")
library("gstat")
library("alexmisc")
library("sp")
library("spgwr")

# Parallelization:
library("future")
library("future.callr")
plan(future.callr::callr, workers = 4)
source("spdiagnostics-functions.R", encoding = "UTF-8")

# Preprocessed data and other objects:
load("run_meuse.rda")

NPERM_LOO <- NPERM_cv


res_lm_loo <- sperrorest(fo, d, model_fun = lm, 
                     coords = c("utmx", "utmy"),
                        smp_fun = partition_loo,
                        smp_args = list(seed1 = 123),
                        importance = TRUE, imp_permutations = NPERM_LOO,
                        imp_sample_from = "all",
                        distance = TRUE, verbose = 1, progress = FALSE,
                        err_fun = err_default)

res_gwr_loo <- sperrorest(fo, d, model_fun = mygwr, 
                      model_args = list(coords = ~utmx+utmy),
                      coords = c("utmx", "utmy"),
                      smp_fun = partition_loo,
                      smp_args = list(seed1 = 123),
                      importance = TRUE, imp_permutations = NPERM_LOO,
                      imp_sample_from = "all",
                      distance = TRUE, verbose = 1, progress = FALSE,
                      mode_rep = "loop", mode_fold = "loop",
                      err_fun = err_default)

res_krg_loo <- sperrorest(fo, d, model_fun = mykrige,
                      model_args = list(range = 500,
                                        locations = ~utmx+utmy),
                      coords = c("utmx", "utmy"),
                      pred_args = list(nmax = NMAX),
                      smp_fun = partition_loo,
                      smp_args = list(seed1 = 123),
                      importance = TRUE, imp_permutations = NPERM_LOO,
                      imp_sample_from = "all",
                      distance = TRUE, verbose = 1, progress = FALSE,
                      mode_rep = "loop", mode_fold = "loop",
                      err_fun = err_default)

res_ok_loo <- sperrorest(fo, d, model_fun = mykrigeOK,
                     model_args = list(range = 500, 
                                       locations = ~ utmx + utmy),
                     pred_args = list(nmax = NMAX),
                     coords = c("utmx", "utmy"),
                     smp_fun = partition_loo,
                     smp_args = list(seed1 = 123),
                     importance = TRUE, imp_permutations = NPERM_LOO,
                     imp_sample_from = "all",
                     distance = TRUE, verbose = 1, progress = FALSE,
                     mode_rep = "loop", mode_fold = "loop",
                     err_fun = err_default)

res_nn_loo <- sperrorest(fo, d, model_fun = nn, 
                     coords = c("utmx", "utmy"),
                     smp_fun = partition_loo,
                     smp_args = list(seed1 = 123),
                     importance = TRUE, imp_permutations = NPERM_LOO,
                     imp_sample_from = "all",
                     distance = TRUE, verbose = 1, progress = FALSE,
                     mode_rep = "loop", mode_fold = "loop",
                     err_fun = err_default)

res_rf_loo <- sperrorest(fo, d, model_fun = randomForest, 
                     model_args = list(ntree = NTREE),
                     coords = c("utmx", "utmy"),
                     smp_fun = partition_loo,
                     smp_args = list(seed1 = 123),
                     importance = TRUE, imp_permutations = NPERM_LOO,
                     imp_sample_from = "all",
                     distance = TRUE, verbose = 1, progress = FALSE,
                     mode_rep = "future", mode_fold = "loop",
                     err_fun = err_default)

res_orf_loo <- sperrorest(fo, d, model_fun = okrf, 
                           model_args = list(ntree = NTREE,
                                             locations = ~utmx+utmy, 
                                             range = 500,
                                             taper_radius = 500),
                           coords = c("utmx", "utmy"),
                           smp_fun = partition_kmeans,
                           smp_args = list(nfold = 10, repetition = 1:NREP_cv, seed1 = 123),
                           importance = TRUE, imp_permutations = NPERM_cv,
                           imp_sample_from = "all",
                           distance = TRUE, verbose = 1, progress = FALSE,
                           mode_rep = "loop", mode_fold = "loop",
                           err_fun = err_default)

save(res_lm_loo, 
     res_gwr_loo,
     res_krg_loo,
     res_ok_loo,
     res_nn_loo,
     res_rf_loo,
     res_orf_loo,
     file = gsub("_res", "_res_loo", OUTFILE))

print(summary(res_rf_loo$error_rep)[c("train_rmse", "test_rmse"), 1:2])
print(summary(res_orf_loo$error_rep)[c("train_rmse", "test_rmse"), 1:2])
print(summary(res_ok_loo$error_rep)[c("train_rmse", "test_rmse"), 1:2])

cat("\nDone.\n")
