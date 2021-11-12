####################################################
### Spatial model diagnostics: Meuse
####################################################
### Resubstitution error
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

# Resampling arguments:
smp_args <- list(nfold = 10, repetition = 1:NREP_cv, seed1 = 123)

# For reproducibility:
set.seed(123)

# List for collecting results:
res_train <- list()

# Coordinates:
loc <- ~ utmx + utmy

# RF:
fit <- randomForest(fo, d, ntree = NTREE)
pred <- predict(fit, newdata = d, nmax = NMAX, locations = loc)
res_train$RF <- err_default(d$logZn, pred)

# OK-RF:
fit <- okrf(fo, d, ntree = NTREE,
            locations = ~utmx+utmy, 
            range = 500,
            taper_radius = 500)
pred <- predict(fit, newdata = d, nmax = NMAX, locations = loc)
res_train$OKRF <- err_default(d$logZn, pred)

# KED:
fit <- mykrige(fo, d, range = 500, locations = loc)
pred <- predict(fit, newdata = d, nmax = NMAX)
res_train$KED <- err_default(d$logZn, pred)

# OK:
fit <- mykrigeOK(fo, d, range = 500, locations = loc)
pred <- predict(fit, newdata = d, nmax = NMAX)
res_train$OK <- err_default(d$logZn, pred)

# LM:
fit <- lm(fo, d)
pred <- predict(fit, newdata = d)
res_train$LM <- err_default(d$logZn, pred)

# GWR:
# Note: Warning related to projection is not relevant.
fit <- mygwr(fo, d, coords = loc)
pred <- predict(fit, newdata = d)
res_train$GWR <- err_default(d$logZn, pred)

# NN:
fit <- nn(fo, d)
pred <- predict(fit, newdata = d)
res_train$NN <- err_default(d$logZn, pred)

# Save results:
save(res_train, file = gsub("_res", "_res_train", OUTFILE))

# Display results:
RES <- res_train %>% dplyr::bind_rows() %>% as.data.frame()
rownames(RES) <- names(res_train)
options(digits = 3)
print(round(RES,3))

cat("\nDone.\n")
