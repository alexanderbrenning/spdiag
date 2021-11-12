####################################################
### Spatial model diagnostics: Maipo
####################################################
### spatial LOO for random forest SPEP
####################################################
### (c) 2021 Alexander Brenning
####################################################

library("sperrorest")
library("randomForest")
library("purrr")

# Additional functions that are not (yet) included in packages:
source("spdiagnostics-functions.R", encoding = "UTF-8")

# Load case study data and associated objects
# saved by spdiag_maipo.R:
load("run_maipo.rda")

res_rf <- sperrorest(fo, d, 
                     model_fun = randomForest,
                     model_args = list(ntree = NTREE),
                     smp_fun = partition_discs4,
                     smp_args = list(maxdist = MAXDIST, 
                                     repetition = 1:NSPLIT, 
                                     seed1 = 123,
                                     n_per_fac = NGROUPS_PER_STRAT,
                                     ndisc = round(NREP/NSPLIT)),
                     distance = FALSE, verbose = 1, progress = FALSE,
                     # do not run this in parallel mode:
                     mode_rep = "loop", mode_fold = "loop",
                     do_gc = 2,
                     err_fun = err_maipo)

save(res_rf, file = gsub("_res", "_res_rf", OUTFILE))
cat("\nDone.\n")

# Not a meaningful summary in this context, just to check that
# the call produced some sort of results:
print(summary(res_rf$error_rep)[c("train_error", "test_error"), c("mean", "sd")])
