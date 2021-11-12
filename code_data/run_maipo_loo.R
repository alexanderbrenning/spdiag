####################################################
### Spatial model diagnostics: Maipo
####################################################
### Leave-one-(grid cell)-out CV
####################################################
### (c) 2021 Alexander Brenning
####################################################

library("sperrorest")
library("randomForest")
library("MASS")
library("purrr")
library("dplyr")
library("magrittr")
library("alexmisc")
library("future")
library("future.callr")

# Parallel processing:
plan(future.callr::callr, workers = 5)

# Additional functions that are not (yet) included in packages:
source("spdiagnostics-functions.R", encoding = "UTF-8")

# Load case study data and associated objects
# saved by spdiag_maipo.R:
load("run_maipo.rda")

# x and y must be passed to splda, will only be used by the
# nearest-neighbour component of NN-LDA, not by the LDA component
spfo <- as.formula(paste(all.vars(fo)[1], "~",
                         as.character(fo)[3], "+ x + y"))

# Arguments for partition_loo_fac_strat:
smp_args <- list(ndisc = round(NREP/NSPLIT), 
                 repetition = 1:NSPLIT, seed1 = 123,
                 strat_fac = "croptype",
                 group_fac = "field",
                 n_per_group = NGROUPS_PER_STRAT)

res_splda_loo <- sperrorest(spfo, d, 
                            model_fun = splda,
                            model_args = list(maxdist = 100),
                            pred_fun = predict.splda,
                            smp_fun = partition_loo_fac_strat,
                            smp_args = smp_args,
                            distance = TRUE, 
                            mode_rep = "loop", mode_fold = "loop",
                            verbose = 1, progress = FALSE)

try(print(summary(res_splda_loo$error_rep)[c("train_error", "test_error"),"mean"]))

res_lda_loo <- sperrorest(fo, d, 
                          model_fun = lda,
                          pred_fun = predict_lda,
                          smp_fun = partition_loo_fac_strat,
                          smp_args = smp_args,
                          distance = TRUE, 
                          mode_rep = "loop", mode_fold = "loop",
                          verbose = 1, progress = FALSE)

try(print(summary(res_lda_loo$error_rep)[c("train_error", "test_error"),"mean"]))

res_rf_loo <- sperrorest(fo, d, 
                         model_fun = randomForest,
                         model_args = list(ntree = NTREE),
                         smp_fun = partition_loo_fac_strat,
                         smp_args = smp_args,
                         distance = TRUE, 
                         mode_rep = "loop", mode_fold = "loop",
                         verbose = 1, progress = FALSE)

try(print(summary(res_rf_loo$error_rep)[c("train_error", "test_error"),"mean"]))

save(res_rf_loo, res_lda_loo, res_splda_loo, 
     file = gsub("_res", "_res_loo", OUTFILE))

cat("\nDone.\n")
