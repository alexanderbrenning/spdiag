####################################################
### Spatial model diagnostics: Maipo
####################################################
### Digest and plot results from Maipo case study
####################################################
### (c) 2021 Alexander Brenning
####################################################

# Load packages:
library("sperrorest")
library("dplyr")
library("magrittr")
library("purrr")
library("alexmisc")
#library("car")
library("wiml")
library("future")
library("future.callr")

# Enable parallel execution of some steps:
plan(future.callr::callr, workers = 5)

# Additional functions that are not (yet) included in packages:
source("spdiagnostics-functions.R", encoding = "UTF-8")


####################################################
### Load data and results files
####################################################

# Load case study data and associated objects
# saved by spdiag_maipo.R:
load("run_maipo.rda")
load("wrp_maipo.rda")

# Load results written by run_maipo_*.R scripts:
load("maipo_res_rf.rda")
load("maipo_res_lda.rda")
load("maipo_res_splda.rda")
load("maipo_res_cv.rda")
load("maipo_res_ncv.rda")
load("maipo_res_spcv.rda")
load("maipo_res_loo.rda")
load("maipo_res_wrf.rda")
load("maipo_res_wlda.rda")
load("maipo_res_wsplda.rda")



####################################################
### Extract CV results (not profiles)
####################################################

### Extract estimates of error rates obtained with different
### types of cross-validation:

mer_lda_cv <- summary(res_lda_cv$error_rep)["test_error", "mean"]
mer_splda_cv <- summary(res_splda_cv$error_rep)["test_error", "mean"]
mer_rf_cv  <- summary(res_rf_cv$error_rep)["test_error", "mean"]

mer_lda_ncv <- summary(res_lda_ncv$error_rep)["test_error", "mean"]
mer_splda_ncv <- summary(res_splda_ncv$error_rep)["test_error", "mean"]
mer_rf_ncv  <- summary(res_rf_ncv$error_rep)["test_error", "mean"]

mer_lda_spcv <- summary(res_lda_spcv$error_rep)["test_error", "mean"]
mer_splda_spcv <- summary(res_splda_spcv$error_rep)["test_error", "mean"]
mer_rf_spcv  <- summary(res_rf_spcv$error_rep)["test_error", "mean"]

mer_lda_loo <- summary(res_lda_loo$error_rep)["test_error", "mean"]
mer_splda_loo <- summary(res_splda_loo$error_rep)["test_error", "mean"]
mer_rf_loo  <- summary(res_rf_loo$error_rep)["test_error", "mean"]

mer_lda_train <- summary(res_lda_cv$error_rep)["train_error", "mean"]
mer_splda_train <- summary(res_splda_cv$error_rep)["train_error", "mean"]
mer_rf_train  <- summary(res_rf_cv$error_rep)["train_error", "mean"]



####################################################
### Extract spatial variable importance profiles
####################################################

# Names of predictors (or transformed features) used in
# permutation-based variable importance assessments:

imp_names <- rownames(summary(res_wlda$importance)) # Early1 etc.


### Spatial variable importance profiles:

# For each of the models:
# - Extract raw importances at the level of individual test points
#   (since test sets are of size 1).
# - Extract associated minimum distances to the respective training sample.
# Result: data frame with one column per (transformed) feature (e.g., 
# Early1), and one additional column 'dist'.

imp_wlda <- res_wlda$importance %>% map(impdists) %>% dplyr::bind_rows()
imp_wlda$dist <- (res_wlda$error_fold %>% map(merdists) %>% dplyr::bind_rows())$dist

imp_wsplda <- res_wsplda$importance %>% map(impdists) %>% dplyr::bind_rows()
imp_wsplda$dist <- (res_wsplda$error_fold %>% map(merdists) %>% dplyr::bind_rows())$dist

imp_wrf <- res_wrf$importance %>% map(impdists) %>% dplyr::bind_rows()
imp_wrf$dist <- (res_wrf$error_fold %>% map(merdists) %>% dplyr::bind_rows())$dist



####################################################
### Extract spatial prediction error profiles
####################################################

# Create data.frames with columns representing error and prediction
# distance. Note that each test fold contains only one observation,
# therefore the error is either 0 or 1.
ed_lda <- res_lda$error_fold %>% map(merdists) %>% dplyr::bind_rows()
ed_splda <- res_splda$error_fold %>% map(merdists) %>% dplyr::bind_rows()
ed_rf <- res_rf$error_fold %>% map(merdists) %>% dplyr::bind_rows()



####################################################
### Create (smoothed) importance and error profiles
####################################################

# Use quadratically increasing prediction distances
# to average the prediction error and variable importance:
brk <- seq(sqrt(30.1), sqrt(MAXDIST), length = 50)^2

# Compute error profiles:

mer_lda <- mers(ed_lda, breaks = brk, smoother = "MA")
mer_splda <- mers(ed_splda, breaks = brk, smoother = "MA")
mer_rf <- mers(ed_rf, breaks = brk, smoother = "MA")

# Compute variable importance profiles:

dmer_wlda <- imp_names %>% 
  map(function(x) mers(imp_wlda[,c(x,"dist")], which = x, breaks = brk, negate = TRUE))
names(dmer_wlda) <- imp_names
dmer_wlda <- dmer_wlda %>% map(2) %>% dplyr::bind_cols() %>% as.data.frame()
### Could use this to rescale importances to add up to 1:
# dmer_wlda <- dmer_wlda / rowSums(dmer_wlda)
dmer_wlda$dist <- mer_lda$dist

dmer_wsplda <- imp_names %>% 
  map(function(x) mers(imp_wsplda[,c(x,"dist")], which = x, breaks = brk, negate = TRUE))
names(dmer_wsplda) <- imp_names
dmer_wsplda <- dmer_wsplda %>% map(2) %>% dplyr::bind_cols() %>% as.data.frame()
# dmer_wsplda <- dmer_wlda / rowSums(dmer_wlda)
dmer_wsplda$dist <- mer_lda$dist

dmer_wrf <- imp_names %>% 
  map(function(x) mers(imp_wrf[,c(x,"dist")], which = x, breaks = brk, negate = TRUE))
names(dmer_wrf) <- imp_names
dmer_wrf <- dmer_wrf %>% map(2) %>% dplyr::bind_cols() %>% as.data.frame()
# dmer_wrf <- dmer_wrf / rowSums(dmer_wrf)
dmer_wrf$dist <- mer_lda$dist


####################################################
### Compile results for plotting in RMarkdown paper
####################################################

# Predictor names:
xvars <- all.vars(fo)[-1]

# Plot parameters:
pal <- list(
  LDA = list(
    name = "LDA",
    col = "blue",
    lty = "solid"
  ),
  SPLDA = list(
    name = "spLDA",
    col = "deepskyblue1",
    lty = "solid"
  ),
  RF = list(
    name = "RF",
    col = "forestgreen",
    lty = "solid"
  )
)

# Variable importance profiles:
dmer <- list(
  SPLDA = dmer_wsplda,
  LDA = dmer_wlda,
  RF = dmer_wrf
)


# Data, results etc.:
MAIPO <- list(
  formula = fo,                    # model formula
  data = d,                        # data from all fields
  xvars = xvars,                   # names of predictors
  nrep = NREP,                     # number of LOO points
  nsplit = NSPLIT,                 # number of times d is resampled
  ntree = NTREE,                   # number of random forest trees
  nfields_per_croptype = NGROUPS_PER_STRAT, # fields sampled per crop type
  maxdist = MAXDIST,
  nn_maxdist = 100,
  mer_lda_cv = mer_lda_cv,         # error rate for LDA in random CV
  mer_lda_loo = mer_lda_loo,       # etc.
  mer_lda_ncv = mer_lda_ncv,
  mer_lda_spcv = mer_lda_spcv,
  mer_lda_train = mer_lda_train,
  mer_lda = mer_lda,
  mer_splda_cv = mer_splda_cv,         
  mer_splda_loo = mer_splda_loo,       
  mer_splda_ncv = mer_splda_ncv,
  mer_splda_spcv = mer_splda_spcv,
  mer_splda_train = mer_splda_train,
  mer_splda = mer_splda,
  mer_rf_cv = mer_rf_cv,
  mer_rf_loo = mer_rf_loo,
  mer_rf_ncv = mer_rf_ncv,
  mer_rf_spcv = mer_rf_spcv,
  mer_rf_train = mer_rf_train,
  mer_rf = mer_rf,
  dmer_wlda = dmer_wlda,           # variable importance profiles for LDA
  dmer_wrf = dmer_wrf,             # variable importance profiles for RF
  dmer = dmer,
  wrp = wrp,                       # transformation object (structured PCA)
  pal = pal,
  NULL
)



####################################################
### Compute minimum distance statistics
### for prediction distance histograms
####################################################

# Distributions of minimum prediction distances
# for different CV settings (for histograms):

getpreddist0 <- function(rsmp) {
  train <- MAIPO$data[rsmp$train, c("x","y") ]
  test <- MAIPO$data[rsmp$test, c("x","y")]
  mindi <- function(i) {
    min(sqrt((train$x - test$x[i])^2 + 
               (train$y - test$y[i])^2))
  }
  preddist <- sapply(1:nrow(test), mindi)
  preddist
}

getpreddist <- function(rsmp) {
  sel <- c(rsmp$train, rsmp$test)
  sel <- c(1:nrow(MAIPO$data)) %in% sel
  train <- MAIPO$data[sel, c("x","y") ]
  test <- MAIPO$data[!sel, c("x","y")]
  mindi <- function(i) {
    min(sqrt((train$x - test$x[i])^2 + 
               (train$y - test$y[i])^2))
  }
  preddist <- sapply(sample(1:nrow(test), 1000), mindi)
  preddist
}

MAIPO$ncvdist <- res_lda_ncv$represampling %>% map(function(x) x[[1]] %>% getpreddist0()) %>% flatten() %>% unlist()
MAIPO$cvdist <- res_lda_cv$represampling %>% map(function(x) x[[1]] %>% getpreddist0()) %>% flatten() %>% unlist()
MAIPO$spcvdist <- res_lda_spcv$represampling %>% map(function(x) x[[1]] %>% getpreddist0()) %>% flatten() %>% unlist()


# Mean minimum predictions distance (one value for each setting):

MAIPO$mean_cvdist <- mean(MAIPO$cvdist)
MAIPO$mean_ncvdist <- mean(MAIPO$ncvdist)
MAIPO$mean_spcvdist <- mean(MAIPO$spcvdist)
MAIPO$mean_loodist <- 30


# Distribution of prediction distances in the prediction scenario:

preddist <- res_lda_cv$represampling %>% 
  map(function(x) x[[1]] %>% getpreddist()) %>% 
  flatten() %>% unlist()
MAIPO$preddist <- preddist
# Summary statistics:
MAIPO$mean_preddist <- mean(preddist)
MAIPO$median_preddist <- median(preddist)
MAIPO$q1_preddist <- quantile(preddist, 0.25)
MAIPO$q3_preddist <- quantile(preddist, 0.75)



####################################################
### Save results for RMarkdown paper
####################################################

saveRDS(MAIPO, file = "maipo_smry.rds")
