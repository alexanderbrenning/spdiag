####################################################
### Spatial model diagnostics: Meuse
####################################################
### Digest results from Meuse case study
####################################################
### (c) 2021 Alexander Brenning
####################################################

# Load packages:
library("sperrorest")
library("gstat")
library("randomForest")
library("sp")
library("purrr")
library("magrittr")
library("alexmisc")
library("future")
library("future.callr")

# Enable parallel execution of some steps:
plan(future.callr::callr, workers = 10)

# Additional functions that are not (yet) included in packages:
source("spdiagnostics-functions.R", encoding = "UTF-8")


####################################################
### Load data and results files
####################################################

# Load case study data and associated objects
# saved by spdiag_meuse.R:
load("run_meuse.rda")

# Load results written by run_meuse_*.R scripts:
# - Distance-based error estimates for each model:
load("meuse_res_rf.rda")
load("meuse_res_nn.rda")
load("meuse_res_lm.rda")
load("meuse_res_krg.rda")
load("meuse_res_ok.rda")
load("meuse_res_gwr.rda")
load("meuse_res_okrf.rda")
# - Resubstitution error:
load("meuse_res_train.rda")
# - CV results using different CV methods:
load("meuse_res_loo.rda")
load("meuse_res_cv.rda")
load("meuse_res_spcv.rda")

# Compile all results in large lists:
res <- list(
  NN = res_nn,
  KED = res_krg,
  OK = res_ok,
  LM = res_lm,
  RF = res_rf,
  OKRF = res_orf,
  GWR = res_gwr
)
rm(res_nn, res_krg, res_ok, res_lm, res_gwr, res_orf)

res_loo <- list(
  NN = res_nn_loo,
  KED = res_krg_loo,
  OK = res_ok_loo,
  LM = res_lm_loo,
  RF = res_rf_loo,
  OKRF = res_orf_loo,
  GWR = res_gwr_loo
)

res_cv <- list(
  NN = res_nn_cv,
  KED = res_krg_cv,
  OK = res_ok_cv,
  LM = res_lm_cv,
  RF = res_rf_cv,
  OKRF = res_orf_cv,
  GWR = res_gwr_cv
)

res_spcv <- list(
  NN = res_nn_spcv,
  KED = res_krg_spcv,
  OK = res_ok_spcv,
  LM = res_lm_spcv,
  RF = res_rf_spcv,
  OKRF = res_orf_spcv,
  GWR = res_gwr_spcv
)


####################################################
### Extract CV results (not profiles)
####################################################

### Extract estimates of RMSE obtained with different
### types of cross-validation:

getCVrmse <- function(x) summary(x$error_rep)["test_rmse","mean"]

rmse_loo <- res_loo %>% map(getCVrmse) %>% unlist()
rmse_cv <- res_cv %>% map(getCVrmse) %>% unlist()
rmse_spcv <- res_spcv %>% map(getCVrmse) %>% unlist()

# Resubstitution RMSE:
rmse_train <- sapply(res_train, function(x) x$rmse)


####################################################
### Extract variable importances (not profiles)
####################################################

getCVimp <- function(obj) {
  res <- summary(obj$importance)
  nms <- rownames(res)
  res <- res[,c("mean.rmse")]
  names(res) <- nms
  res
}

imp_loo <- res_loo %>% map(getCVimp) %>% dplyr::bind_rows() %>% as.data.frame()
rownames(imp_loo) <- names(res_loo)
imp_cv <- res_cv %>% map(getCVimp) %>% dplyr::bind_rows() %>% as.data.frame()
rownames(imp_cv) <- names(res_cv)
imp_spcv <- res_spcv %>% map(getCVimp) %>% dplyr::bind_rows() %>% as.data.frame()
rownames(imp_spcv) <- names(res_spcv)


####################################################
### Extract spatial prediction error profiles
####################################################

extract_errors <- function(res, type = "test") {
  res$error_fold %>% map(errdists, type = type) %>% dplyr::bind_rows()
}

# Spatial leave-one-out RMSE:
ed <- lapply(res, extract_errors)

# Compute smoothed prediction error profiles:
brk <- seq(20, MAXDIST, length = 30)
rmse <- lapply(ed, function(x) smth(rmses0(x, breaks = brk)))



####################################################
### Extract variable importances profiles
####################################################

# This is a bit complex because sperrorest currently only returns
# differences between permuted performance estimated and unpermuted
# performance estimates. In LOO estimation, there is only one
# observation in each test set, and therefore the original permuted
# predictions can be reconstructed. These are afterwards used
# to calculate RMSEs from permuted and unpermuted predictions.

# Spatial importance profiles:

extract_importance <- function(res) {
  # In this context, 'bias' is simply the difference between observed
  # and predicted - there is no averaging since the test set contains
  # only one observation in LOO CV.
  imp <- res$importance[ sapply(res$importance, function(x) class(x) != "resampling") ]
  imp <- imp %>% flatten() %>% map("bias") %>% as.data.frame() %>% t()
  rownames(imp) <- NULL
  vnms <- rownames(res$importance[[1]][[1]])
  colnames(imp) <- vnms
  imp <- as.data.frame(imp)
  imp$dist <- (res$error_fold %>% map(biasdists) %>% dplyr::bind_rows())$dist
  imp$obs <- res$error_fold %>% map(function(x) x[[1]]$test$obs) %>% unlist() %>% unname()
  imp$pred <- res$error_fold %>% map(function(x) x[[1]]$test$pred) %>% unlist() %>% unname()
  imp$bias <- res$error_fold %>% map(function(x) x[[1]]$test$bias) %>% unlist() %>% unname()
  # Reconstruct the value predicted with permuted feature values:
  for (vnm in vnms) {
    imp[,vnm] <- imp[,vnm] + imp$bias
  }
  imp
}

extract_dmer <- function(imp, breaks, nms = NULL, normalize = FALSE) {
  if (is.null(nms)) {
    nms <- colnames(imp)
    nms <- nms[ !(nms %in% c("dist", "obs", "pred", "bias")) ]
  }
  dmer <- nms %>% 
    map(function(x) rmses0(imp[,c(x,"obs","dist")], which = x, breaks = breaks))
  names(dmer) <- nms
  dmer <- dmer %>% map(2) %>% dplyr::bind_cols() %>% as.data.frame()
  if (normalize)
    dmer <- dmer / rowSums(dmer)
  dmer$dist <- (breaks[-1] + breaks[-length(breaks)]) / 2
  dmer <- smth(dmer)
  dmer
}

# For each model, extract permuted predictions and distances:
imp <- lapply(res, extract_importance)
# These are NOT yet the importances.

# For each model, compute binned permuted predictions:
dmer <- lapply(imp, extract_dmer, breaks = brk)

# Compute permutation importances by subtracting the corresponding
# unpermuted RMSE:
vnms <- all.vars(fo)[-1]
for (m in names(dmer)) { # for each model
  dmer[[m]] <- dmer[[m]][,c(vnms, "dist")] # drop 'dist'
  for (vnm in vnms) { # for each predictor
    dmer[[m]][,vnm] <- dmer[[m]][,vnm] - rmse[[m]]$rmse
  }
  # Smoothed profile:
  dmer[[m]] <- smth(dmer[[m]])
}


####################################################
### Compute training-set sample sizes
####################################################

# This is just a side note... need to know how rapidly the training
# sample sizes decrease as the buffer distance increases.
# Result is briefly mentioned in the paper.

extract_counts <- function(x) data.frame(train_count = x[[1]]$train$count, test_count = x[[1]]$test$count, dist = x[[1]]$distance)

count <- res$NN$error_fold %>% map(extract_counts) %>% dplyr::bind_rows() %>% as.data.frame()



####################################################
### Extract variable importance profiles
####################################################

# Compile data and results for use in the RMarkdown
# document that generates the paper:

MEUSE <- list(
  data = d,        # data
  formula = fo,    # model formula
  nrep = NREP,     # number of repetitions, i.e. LOO points
  nperm = NPERM,   # number of permutations (=1)
  nmax = NMAX,     # number of nearest neighbours in IDW, OK, KED
  ntree = NTREE,   # number of bagged trees
  rmse = rmse,     # smoothed RMSE profiles for plotting
  dmer = dmer,     # smoothed SVI profiles for plotting
  rmse_train = rmse_train,  # resubstitution error
  rmse_loo = rmse_loo,
  imp_loo = imp_loo,
  rmse_spcv = rmse_spcv,
  imp_spcv = imp_spcv,
  rmse_cv = rmse_cv,
  imp_cv = imp_cv,
  count = count    # sample sizes from each LOO run
)


####################################################
### Compute prediction distances
####################################################

# Minimum nearest-neighbor distance:
di <- as.matrix(dist(MEUSE$data[,c("utmx","utmy")]))
diag(di) <- NA
MEUSE$mindist <- apply(di, 1, min, na.rm=TRUE)
MEUSE$mean_mindist <- mean(MEUSE$mindist)
rm(di)


getpreddist0 <- function(rsmp) {
  train <- MEUSE$data[rsmp$train, c("utmx","utmy") ]
  test <- MEUSE$data[rsmp$test, c("utmx","utmy")]
  mindi <- function(i) {
    min(sqrt((train$utmx - test$utmx[i])^2 + 
             (train$utmy - test$utmy[i])^2))
  }
  preddist <- sapply(1:nrow(test), mindi)
  preddist
}

getpreddist <- function(x) {
  x$represampling %>% 
  map(function(y) y %>% map(getpreddist0) %>% flatten()) %>% 
  flatten() %>% unlist() %>% unname()
}

# Prediction distances in various types of cross-validation
# for histograms in the Appendix:
MEUSE$loodist <- getpreddist(res_lm_loo)
MEUSE$cvdist <- getpreddist(res_lm_cv)
MEUSE$spcvdist <- getpreddist(res_lm_spcv)

# Mean prediction distances in various types of CV:
MEUSE$mean_cvdist <- mean(MEUSE$cvdist)
MEUSE$mean_loodist <- mean(MEUSE$loodist)
MEUSE$mean_spcvdist <- mean(MEUSE$spcvdist)

# Prediction distances when prediction on the entire
# Meuse floodplain:
data("meuse.grid", package = "sp")
mindi <- function(i) {
  min(sqrt((MEUSE$data$utmx - meuse.grid$x[i])^2 + 
           (MEUSE$data$utmy - meuse.grid$y[i])^2))
}
preddist <- sapply(seq_len(nrow(meuse.grid)), mindi)
MEUSE$preddist <- preddist
MEUSE$mean_preddist <- mean(preddist)
MEUSE$median_preddist <- median(preddist)
MEUSE$q1_preddist <- quantile(preddist, 0.25)
MEUSE$q3_preddist <- quantile(preddist, 0.75)


####################################################
### Regression and semivariogram analysis
####################################################

# Just some preliminary analyses to inform the reader
# about the strength of the predictor-response relationship
# and the spatial structure of the response and residuals.

# Multiple linear regression:
fit <- lm(MEUSE$formula, MEUSE$data)
MEUSE$r.squared <- summary(fit)$r.squared

# Semivariogram analysis:
library("sp")
library("gstat")
meu <- MEUSE$data
coordinates(meu) <- ~utmx+utmy
vg <- variogram(logZn~1, meu)
fit_vg <- fit.variogram(vg, vgm(1, "Sph", 500, 1))
rvg <- variogram(MEUSE$formula, meu)
fit_rvg <- fit.variogram(rvg, vgm(1, "Sph", 500, 1))
MEUSE$svgm <- list(
  nugget = fit_vg$psill[1],
  sill = sum(fit_vg$psill),
  ns_ratio = fit_vg$psill[1] / sum(fit_vg$psill),
  range = fit_vg$range[2]
)
MEUSE$resid_svgm <- list(
  nugget = fit_rvg$psill[1],
  sill = sum(fit_rvg$psill),
  ns_ratio = fit_rvg$psill[1] / sum(fit_rvg$psill),
  range = fit_rvg$range[2]
)
rm(fit, meu, vg, rvg, fit_vg, fit_rvg)


####################################################
### Prediction maps of selected models for Figure
####################################################

library("sp")
library("gstat")
library("raster")

# Prepare data.frame for prediciton on grid:
MEUSE$newdata <- meuse.grid
MEUSE$newdata$utmx <- MEUSE$newdata$x
MEUSE$newdata$utmy <- MEUSE$newdata$y
MEUSE$newdata$x <- MEUSE$newdata$y <- NULL
MEUSE$newdata$sqrt.dist <- sqrt(MEUSE$newdata$dist)

# Interpolate digital elevation model:
elev_ked <- mykrige(elev~sqrt.dist, d, locations = ~utmx+utmy, range = 500)
MEUSE$newdata$elev <- predict(elev_ked, newdata = MEUSE$newdata)

# Random forest prediction:
fit <- randomForest(fo, MEUSE$data, ntree = MEUSE$ntree)
MEUSE$newdata$pred_rf <- predict(fit, newdata = MEUSE$newdata)

# Kriging with external drift:
fit <- mykrige(fo, MEUSE$data, locations = ~utmx+utmy, range = 500)
MEUSE$newdata$pred_ked <- predict(fit, newdata = MEUSE$newdata)

# Ordinary kriging:
fit <- mykrigeOK(fo, MEUSE$data, locations = ~utmx+utmy, range = 500)
MEUSE$newdata$pred_ok <- predict(fit, newdata = MEUSE$newdata)

# Hybrid OK-RF method:
fit <- okrf(fo, MEUSE$data, locations = ~utmx+utmy, range = 500)
MEUSE$newdata$pred_okrf <- predict(fit, newdata = MEUSE$newdata)

# Convert predictions to raster objects:
MEUSE$spnewdata <- MEUSE$newdata
coordinates(MEUSE$spnewdata) <- ~utmx+utmy
MEUSE$spnewdata <- as(MEUSE$spnewdata, "SpatialPixelsDataFrame")
MEUSE$r_pred_rf <- raster::rasterize(MEUSE$spnewdata, 
                                     raster::raster(as(MEUSE$spnewdata, "SpatialPixels")), 
                                     field = "pred_rf")
MEUSE$r_pred_ked <- raster::rasterize(MEUSE$spnewdata, 
                                      raster::raster(as(MEUSE$spnewdata, "SpatialPixels")), 
                                      field = "pred_ked")
MEUSE$r_pred_ok <- raster::rasterize(MEUSE$spnewdata, 
                                     raster::raster(as(MEUSE$spnewdata, "SpatialPixels")), 
                                     field = "pred_ok")
MEUSE$r_pred_okrf <- raster::rasterize(MEUSE$spnewdata, 
                                       raster::raster(as(MEUSE$spnewdata, "SpatialPixels")), 
                                       field = "pred_okrf")




####################################################
### Plot settings
####################################################

pal <- list(
  NN = list(
    name = "NN",
    col = "brown",
    lty = "solid"),
  KED = list(
    name = "KED",
    col = "black",
    lty = "solid"),
  OK = list(
    name = "OK",
    col = "grey60",
    lty = "solid"),
  LM = list(
    name = "LM",
    col = "blue",
    lty = "solid"),
  GWR = list(
    name = "GWR",
    col = "deepskyblue1",
    lty = "solid"),
  OKRF = list(
    name = "OKRF",
    col = "green",
    lty = "solid"),
  RF = list(
    name = "RF",
    col = "forestgreen",
    lty = "solid")
)

MEUSE$pal <- pal



####################################################
### Save results and data for RMarkdown paper
####################################################

saveRDS(MEUSE, file = "meuse_smry.rds")
