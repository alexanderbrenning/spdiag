####################################################
### Spatial model diagnostics: Maipo
####################################################
### spatial LOO for NN-LDA SPEP
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

spfo <- as.formula(paste(all.vars(fo)[1], "~",
                 as.character(fo)[3], "+ x + y"))

res_splda <- sperrorest(spfo, d, model_fun = splda,
                        model_args = list(maxdist = 100),
                        pred_fun = predict.splda,
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

save(res_splda, file = gsub("_res", "_res_splda", OUTFILE))
cat("\nDone.\n")

print(summary(res_splda$error_rep)[c("train_error", "test_error"), c("mean", "sd")])
