####################################################
### Spatial model diagnostics: Maipo
####################################################
### SVIPs of NN-LDA in transformed feature space
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

# x and y must be passed to splda, will only be used by the
# nearest-neighbour component of NN-LDA, not by the LDA component
spfo <- as.formula(paste(all.vars(fo)[1], "~",
                         as.character(fo)[3], "+ x + y"))

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

warped_splda <- warp(splda, warper = wrp)

wspfo <- warp(spfo, warper = wrp)

res_wsplda <- sperrorest(wspfo, wd, 
                         model_fun = warped_splda,
                         model_args = list(maxdist = 100),
                         pred_fun = predict.wsplda,
                         smp_fun = partition_discs4,
                         smp_args = smp_args,
                         distance = FALSE, 
                         verbose = 1, progress = FALSE,
                         err_fun = err_maipo,
                         mode_rep = "loop", mode_fold = "loop",
                         importance = TRUE, imp_permutations = 10,
                         imp_sample_from = "all",
                         imp_variables = imp_variables)

save(res_wsplda, file = gsub("_res", "_res_wsplda", OUTFILE))

print(summary(res_wsplda$error_rep)[c("train_error", "test_error"), c("mean", "sd")])
print(summary(res_wsplda$importance)[, c("mean.error", "sd.error")])

cat("\nDone.\n")
