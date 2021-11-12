####################################################
### Spatial model diagnostics: Maipo
####################################################
### spatial LOO for LDA SPEP
####################################################
### (c) 2021 Alexander Brenning
####################################################

### For detailed comments see run_maipo_rf.R!

library("sperrorest")
library("MASS")
library("purrr")
library("dplyr")
library("magrittr")
library("alexmisc")

source("spdiagnostics-functions.R", encoding = "UTF-8")

load("run_maipo.rda")

# Predict function for LDA that predicts class membership:
predict_lda <- function(object, newdata)predict(object, newdata)$class


res_lda <- sperrorest(fo, d, model_fun = lda,
                      pred_fun = predict_lda,
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

save(res_lda, file = gsub("_res", "_res_lda", OUTFILE))
cat("\nDone.\n")

print(summary(res_lda$error_rep)[c("train_error", "test_error"), c("mean", "sd")])
