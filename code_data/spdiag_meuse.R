####################################################
### Spatial model diagnostics: Meuse
####################################################
### Prepare Meuse data and model settings
####################################################
### (c) 2021 Alexander Brenning
####################################################

library("sp")
data("meuse", package = "sp")
d <- meuse
d$logZn <- log10(d$zinc)
d$sqrt.dist <- sqrt(d$dist)
d$flood <- as.numeric(as.character(d$ffreq))
d$flood <- as.numeric(d$flood > 1)
d$utmx <- d$x
d$utmy <- d$y
d$x <- d$y <- d$dist.m <- NULL

fo <- logZn ~ sqrt.dist + elev + utmx + utmy


MAXDIST <- 1100    # Maximum prediction distance for spatial LOO
NREP <- 50000      # Number of repetitions in spatial LOO
NTREE <- 500       # 
NMAX <- 500 # Inf
NPERM <- 1
NPERM_cv <- 10
NREP_cv <- 10

OUTFILE <- "meuse_res.rda"

save(d, fo, 
     MAXDIST, OUTFILE, NTREE, NREP, NPERM, NMAX,
     NREP_cv, NPERM_cv,
     file = "run_meuse.rda")
