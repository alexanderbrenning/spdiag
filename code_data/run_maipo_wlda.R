####################################################
### Spatial model diagnostics: Maipo
####################################################
### SVIPs of LDA in transformed feature space
####################################################
### (c) 2021 Alexander Brenning
####################################################

library("sperrorest")
library("MASS")
library("purrr")
library("dplyr")
library("magrittr")
library("alexmisc")
library("future")
library("future.callr")
library("iml")
library("wiml")

# Additional functions that are not (yet) included in packages:
source("spdiagnostics-functions.R", encoding = "UTF-8")

# Load case study data and associated objects
# saved by spdiag_maipo.R:
load("run_maipo.rda")
load("wrp_maipo.rda")

# Arguments for resampling function:
smp_args <- list(maxdist = MAXDIST, 
                 repetition = 1:NSPLIT, 
                 seed1 = 123,
                 n_per_fac = NGROUPS_PER_STRAT,
                 ndisc = round(NREP/NSPLIT))

imp_variables <- c("Early1", "Early2", "Early3", 
                   "Early4", "Early5", "Early6",
                   "Mid1", "Mid2", "Mid3", "Mid4",
                   "Late1", "Late2", "Late3", 
                   "Late4", "Late5", "Late6", 
                   "Late7", "Late8")

res_wlda <- sperrorest(wfo, wd, model_fun = warped_lda,
                       pred_fun = predict_lda,
                       smp_fun = partition_discs4,
                       smp_args = smp_args,
                       distance = FALSE, 
                       verbose = 1, progress = FALSE,
                       err_fun = err_maipo,
                       # do not run this in parallel mode:
                       mode_rep = "loop", mode_fold = "loop",
                       importance = TRUE, imp_permutations = 10,
                       imp_sample_from = "all",
                       imp_variables = imp_variables)

save(res_wlda, file = gsub("_res", "_res_wlda", OUTFILE))

print(summary(res_wlda$error_rep)[c("train_error", "test_error"), c("mean", "sd")])
print(summary(res_wlda$importance)[, c("mean.error", "sd.error")])

cat("\nDone.\n")
