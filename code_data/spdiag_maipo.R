####################################################
### Spatial model diagnostics: Maipo
####################################################
### Prepare Maipo data and model settings,
### save to run_Maipo.rda and wrp_Maipo.rda
####################################################
### (c) 2021 Alexander Brenning
####################################################

# Required packages:
library("purrr") # flatten()
library("MASS") # lda()
library("randomForest") # randomForest()
library("wiml") # warp()

# Additional functions such as splda():
source("spdiagnostics-functions.R", encoding = "UTF-8")

data("maipo", package = "sperrorest")
d <- maipo

# Prepare data.frame:
d$field <- droplevels(d$field)
d$x <- d$utmx
d$y <- d$utmy
# Names of predictors and response:
xvars <- stringr::str_subset(colnames(d), "nd.i0.")
xvars <- c(xvars, stringr::str_subset(colnames(d), "^b[0-9]{2}"))
yvar <- "croptype"

# Reduce dataset to required variables:
d <- d[, c(yvar, xvars, "x", "y", "field")]

# Model formula with all predictors:
fo <- as.formula(paste(yvar, "~", paste(xvars, collapse = "+")))

# Cross-validation settings:
MAXDIST <- 15000
NSPLIT <- 50
NREP <- 25000
NGROUPS_PER_STRAT <- N_PER_FIELD <- 25

# Random forest model setting:
NTREE <- 500

# Base file name for results files:
OUTFILE <- "maipo_res.rda"


save(d, fo, MAXDIST, OUTFILE, NTREE, NREP, NSPLIT, 
     N_PER_FIELD, NGROUPS_PER_STRAT,
     file = "run_maipo.rda")


####################################################
### Prepare transformed Maipo data and models
####################################################

yvar <- "croptype"

# Structured list of predictors for structured PCA:
xvars_list <- list(
  "Early" = c("ndvi01", "ndvi02", "ndwi01", "ndwi02", stringr::str_subset(xvars, "^b[1-2][0-9]$")),
  "Mid" = c("ndvi03", "ndwi03", stringr::str_subset(xvars, "^b3[0-9]$")),
  "Late" = c("ndvi04", "ndvi05", "ndvi06", "ndvi07", "ndvi08", 
             "ndwi04", "ndwi05", "ndwi06", "ndwi07", "ndwi08",
             stringr::str_subset(xvars, "^b[4-8][0-9]$")))

# Keep x/y variables as untransformed variables:
uvars_list <- c("x", "y")

# Strucutred PCA transformation object:
wrp <- strucpca_warper(d, xvars = xvars_list, yvar = yvar, 
                       uvars = uvars_list, 
                       wvars = names(xvars_list))

# Apply this transformation to models:
warped_lda <- warp(lda, warper = wrp)
warped_splda <- warp(splda, warper = wrp)
warped_rf <- warp(randomForest, warper = wrp)

# ...and data:
wd <- warp(d, wrp)

### ...and model formula:
wfo <- warp(fo, wrp)

save(yvar, xvars_list, uvars_list, wrp, 
     warped_lda, warped_splda, warped_rf, 
     wd, wfo, 
     file = "wrp_maipo.rda")
