####################################################
### Spatial model diagnostics: Maipo
####################################################
### Spatial CV (k-means coordinate clustering)
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

# Arguments for partition_kmeans_fac_strat:
smp_args <- list(nfold = 10, 
                 repetition = 1:50,
                 strat_fac = "croptype",
                 group_fac = "field",
                 ngroups_per_strat = N_PER_FIELD,
                 balancing_steps = 5,
                 seed1 = 123)

res_splda_spcv <- sperrorest(spfo, d, 
                             model_fun = splda,
                             model_args = list(maxdist = 100),
                             pred_fun = predict.splda,
                             smp_fun = partition_kmeans_fac_strat,
                             smp_args = smp_args,
                             distance = FALSE, 
                             verbose = 1, progress = FALSE)

res_rf_spcv <- sperrorest(fo, d, 
                          model_fun = randomForest,
                          model_args = list(ntree = NTREE),
                          smp_fun = partition_kmeans_fac_strat,
                          smp_args = smp_args,
                          distance = FALSE, 
                          verbose = 1, progress = FALSE)

res_lda_spcv <- sperrorest(fo, d, 
                           model_fun = lda,
                           pred_fun = predict_lda,
                           smp_fun = partition_kmeans_fac_strat,
                           smp_args = smp_args,
                           distance = FALSE, 
                           verbose = 1, progress = FALSE)

save(res_rf_spcv, res_lda_spcv, res_splda_spcv, 
     file = gsub("_res", "_res_spcv", OUTFILE))
cat("\nDone.\n")

print(summary(res_rf_spcv$error_rep)[c("train_error", "test_error"), c("mean", "sd")])
print(summary(res_lda_spcv$error_rep)[c("train_error", "test_error"), c("mean", "sd")])
print(summary(res_splda_spcv$error_rep)[c("train_error", "test_error"), c("mean", "sd")])
