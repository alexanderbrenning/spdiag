####################################################
### Spatial model diagnostics: Maipo
####################################################
### Field-level CV
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

# Arguments for partition_factor_cv_fac_strat:
smp_args <- list(nfold = 10, 
                 repetition = 1:50,
                 strat_fac = "croptype",
                 group_fac = "field",
                 ngroups_per_strat = NGROUPS_PER_STRAT,
                 seed1 = 123)

# Field-level CV of NN-LDA:
res_splda_cv <- sperrorest(spfo, d, 
                           model_fun = splda,
                           model_args = list(maxdist = 100),
                           pred_fun = predict.splda,
                           smp_fun = partition_factor_cv_fac_strat,
                           smp_args = smp_args,
                           mode_rep = "loop", mode_fold = "loop",
                           distance = FALSE, 
                           verbose = 1, progress = FALSE)

# Field-level CV of LDA:
res_lda_cv <- sperrorest(fo, d, 
                         model_fun = lda,
                         pred_fun = predict_lda,
                         smp_fun = partition_factor_cv_fac_strat,
                         smp_args = smp_args,
                         mode_rep = "loop", mode_fold = "loop",
                         distance = FALSE, 
                         verbose = 1, progress = FALSE)

# Field-level CV of RF:
res_rf_cv <- sperrorest(fo, d, 
                        model_fun = randomForest,
                        model_args = list(ntree = NTREE),
                        smp_fun = partition_factor_cv_fac_strat,
                        smp_args = smp_args,
                        mode_rep = "loop", mode_fold = "loop",
                        distance = FALSE, 
                        verbose = 1, progress = FALSE)

save(res_rf_cv, res_lda_cv, res_splda_cv, 
     file = gsub("_res", "_res_cv", OUTFILE))
cat("\nDone.\n")

# Quick summary just to see if it worked:
print(summary(res_rf_cv$error_rep)[c("train_error", "test_error"),"mean"])
print(summary(res_lda_cv$error_rep)[c("train_error", "test_error"),"mean"])
print(summary(res_splda_cv$error_rep)[c("train_error", "test_error"),"mean"])
