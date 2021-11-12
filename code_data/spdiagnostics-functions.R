####################################################
### Spatial model diagnostics: Functions
####################################################
### The functions included in this file are
### mainly from one of the following categories:
### - Partitioning / resampling methods to be used
###   with sperroest
### - Functions for extracting error or variable
###   importance information from sperrorest
###   results objects
### - Error (or loss) functions for sperrorest
### - Customized models such as NN-LDA and OK-RF
####################################################
### (c) 2021 Alexander Brenning
####################################################


####################################################
### Partitioning functions
####################################################


partition_discs <- function(data, coords = c("x", "y"),
                            method = c("interval", "histogram")[2],
                            maxndist = 10,
                            maxdist = Inf,
                            zerodist = 1e-10,
                            repetition = nrow(data),
                            seed1 = NULL) {
  if (length(repetition) < 2) {
    repetition <- seq(1, repetition, by = 1)
  }
  maxndist <- max(2, maxndist)
  maxdist <- max(0, maxdist)
  
  point_discs <- function(cnt, data, coords, method, maxdist, maxndist, seed1) {
    if (!is.null(seed1)) set.seed(seed1 + cnt)
    i <- (cnt %% nrow(data))
    if (i == 0) i <- nrow(data)
    #print(i)
    di <- sqrt((data[, coords[1]] - data[i, coords[1]])^2 + 
                 (data[, coords[2]] - data[i, coords[2]])^2)
    
    if (method == "interval") {
      if (sum(di <= zerodist) >= 2) {
        dists <- c(-1, seq(0, maxdist, length = maxndist - 1))
      } else {
        dists <- seq(0, maxdist, length = maxndist)
      }
      sel <- c(rep(TRUE, length(dists)), FALSE)
      dists <- c(dists, Inf)
      for (i in 2:(length(dists) - 1)) {
        if (!any((di >= dists[i]) & (di < dists[i+1])))
          sel[i] <- FALSE
      }
      dists <- dists[sel]
    } else if (method == "histogram") {
      dists <- sort(unique(di[-i][ di[-i] < maxdist ], decreasing = FALSE))
      #print(dists)
      if (length(dists) == 0) {
        dists <- 0
      } else if (dists[1] <= zerodist) { # if there are other observations with identical coordinate
        dists <- c(-1, dists)
      }
      if (length(dists) > maxndist) {
        dists <- dists[round(seq(1, length(dists), length = maxndist))]
        # dists <- c(dists[1], sample(dists[-1], size = maxndist))
        # dists <- sort(dists, decreasing = FALSE)
      }
    }
    
    repet <- list()
    for (fld in 1:length(dists)) {
      train <- which(di >= dists[fld])
      train <- train[ train != i ]
      repet[[fld]] <- list(train = train, test = c(i))
    }
    repet
  }
  
  resample <- future.apply::future_lapply(repetition, point_discs,
                                          data = data, 
                                          coords = coords, 
                                          method = method,
                                          maxdist = maxdist, 
                                          maxndist = maxndist, 
                                          seed1 = seed1, 
                                          future.seed = TRUE) 
  names(resample) <- as.character(repetition) 
  
  # resample <- list()
  # for (cnt in repetition) {
  #   resample[[as.character(cnt)]] <- repet
  # }
  as.represampling(resample)
}


partition_discs2 <- function(data, coords = c("x", "y"),
                            maxdist = NULL,
                            repetition = min(c(1000, nrow(data)*maxndist)),
                            train_fun = NULL,
                            train_param = NULL,
                            seed1 = NULL) {
  if (length(repetition) < 2) {
    repetition <- seq(1, repetition, by = 1)
  }
  if (is.null(maxdist)) {
    maxdist <- min(c(diff(range(data[,coords[1]])),
                     diff(range(data[,coords[2]])))) * 0.25
  }
  maxdist <- max(0, maxdist)
  
  point_disc <- function(cnt, data, coords, maxdist, 
                         train_fun, train_param, seed1) {
    if (!is.null(seed1)) set.seed(seed1 + cnt)
    if (!is.null(train_fun))
      data <- train_fun(data, param = train_param)
    i <- sample(1:nrow(data), size = 1)
    di <- sqrt((data[, coords[1]] - data[i, coords[1]])^2 +
                 (data[, coords[2]] - data[i, coords[2]])^2)
    # thr <- sample(c(0:maxndist), size = 1) * maxdist / maxndist
    maxdist2 <- maxdist * ifelse(runif(1) < 0.3, 0.3, 1)
    thr <- runif(n = 1, min = 0, max = maxdist2)
    attempts <- 0
    while (all(di > thr) | (all(di <= thr))) {
      maxdist2 <- maxdist * ifelse(runif(1) < 0.3, 0.3, 1)
      thr <- runif(n = 1, min = 0, max = maxdist2)
      attempts <- attempts + 1
      if (attempts >= 20) return(NULL)
    }
    train <- which(di > thr)
    mindi <- min(di[train])
    repet <- list(list(train = train, test = c(i),
                       distance = mindi))
    repet
  }
  
  resample <- lapply(repetition, point_disc,
                     data = data, 
                     coords = coords, 
                     maxdist = maxdist, 
                     seed1 = seed1, 
                     train_fun = train_fun,
                     train_param = train_param) 
  names(resample) <- as.character(repetition) 
  
  as.represampling(resample)
}



partition_discs3 <- function(data, coords = c("x", "y"),
                             maxdist = NULL,
                             ndisc = min(c(200, nrow(data))),
                             repetition = 1,
                             seed1 = NULL,
                             verbose = 0) {
  if (length(repetition) < 2) {
    repetition <- seq(1, repetition, by = 1)
  }
  if (is.null(maxdist)) {
    maxdist <- min(c(diff(range(data[,coords[1]])),
                     diff(range(data[,coords[2]])))) * 0.25
  }
  maxdist <- max(0, maxdist)
  
  if (is.null(seed1)) 
    seed1 <- round(runif(1, 1, 10^9))

  rsmp <- future.apply::future_lapply(as.list(repetition),
                 function(x) {
                   flatten(unclass(
                     partition_discs2(
                       data = data,
                       coords = coords,
                       maxdist = maxdist,
                       repetition = 1:ndisc,
                       seed1 = seed1 + x))) },
                 future.seed = TRUE)
  
  names(rsmp) <- as.character(repetition)

  as.represampling(rsmp)
}


#parti <- partition_discs4(d, ndisc = 10, repetition = 1:3, seed1 = 123)

partition_discs4 <- function(data, coords = c("x", "y"),
                             maxdist = NULL,
                             ndisc = min(c(200, nrow(data))),
                             strat = "croptype",
                             fac = "field",
                             n_per_fac = 25,
                             repetition = 1,
                             seed1 = NULL,
                             verbose = 0) {
  if (length(repetition) < 2) {
    repetition <- seq(1, repetition, by = 1)
  }
  if (is.null(maxdist)) {
    maxdist <- min(c(diff(range(data[,coords[1]])),
                     diff(range(data[,coords[2]])))) * 0.25
  }
  maxdist <- max(0, maxdist)
  
  if (is.null(seed1)) 
    seed1 <- round(runif(1, 1, 10^9))
  
  data$.fid <- seq_len(nrow(data))
  
  rsmp <- list()
  
  for (i in repetition) {
    if (!is.null(seed1))
      set.seed(seed1 + i)
    wh <- c()
    for (lev in levels(data[, strat])) {
      sel_lev <- data[, strat] == lev
      all_fields <- unique(as.character(data[sel_lev, fac]))
      sel_fields <- sample(all_fields, size = n_per_fac, replace = FALSE)
      wh <- c(wh, which(data[, fac] %in% sel_fields))
    }
    dat <- data[wh, ]
    rsmp[[i]] <- flatten(unclass(
      partition_discs2(
        data = dat,
        coords = coords,
        maxdist = maxdist,
        repetition = 1:ndisc)))
    for (j in 1:length(rsmp[[i]])) {
      rsmp[[i]][[j]]$train <- dat$.fid[rsmp[[i]][[j]]$train]
      rsmp[[i]][[j]]$test <- dat$.fid[rsmp[[i]][[j]]$test]
    }
  }

  names(rsmp) <- as.character(repetition)
  
  as.represampling(rsmp)
}




# Helper function for resampling fields in the Maipo case study:
# - Sample at the field level, not at the level of grid cells,
#   which are the actual observations.
# - Retrieve the desired number of fields per crop type, e.g. 50 each
resample_fac_strat <- function(data, param = list(strat_fac = "class", 
                                                  group_fac = "field",
                                                  ngroups_per_strat,
                                                  relevel = TRUE)) {
  if (is.null(param$strat_fac)) {
    param$strat_fac <- "class"
  }
  if (is.null(param$group_fac)) {
    param$group_fac <- "field"
  }
  if (is.null(param$relevel)) {
    param$relevel <- TRUE
  }
  stopifnot(is.numeric(param$ngroups_per_strat))
  
  wh <- c()
  for (lev in levels(data[, param$strat_fac])) {
    sel_lev <- data[, param$strat_fac] == lev
    all_groups <- unique(as.character(data[sel_lev, param$group_fac]))
    sel_groups <- sample(all_groups, size = param$ngroups_per_strat, replace = FALSE)
    wh <- c(wh, which(data[, param$group_fac] %in% sel_groups))
  }
  
  data <- data[wh, ]
  
  if (param$relevel) {
    data[, param$group_fac] <- factor(data[, param$group_fac])
  }
  
  data
}

#dd <- resample_fac_strat(d, param = list(strat_fac = "croptype", group_fac = "field", n_per_group = 25))

# LOO partitioning method for the Maipo case study:
# - A certain number of fields is selected from each crop type.
# - One grid cell at a time is used as the test sample.
# This is non-spatial LOO with a twist to obtain the desired
# training sample size.
partition_loo_fac_strat <- function (data, 
                                     repetition = 1:1,
                                     ndisc = nrow(data), replace = FALSE,
                                     strat_fac = "class", 
                                     group_fac = "field",
                                     n_per_group,
                                     relevel = TRUE,
                                     seed1 = NULL, ...) 
{
  if (length(repetition) < 2) {
    repetition <- seq(1, repetition, by = 1)
  }
  data$.fid <- seq_len(nrow(data))
  rsmp <- list()
  for (i in repetition) {
    if (!is.null(seed1))
      set.seed(seed1 + i)
    dat <- resample_fac_strat(data = data, 
                              param = list(strat_fac = strat_fac, group_fac = group_fac, 
                                           ngroups_per_strat = n_per_group, relevel = relevel))
    rsmp[[i]] <- list()
    for (j in 1:ndisc) {
      loo <- sample(1:nrow(dat), size = 1)
      rsmp[[i]][[j]] <- list(train = dat$.fid[-loo],
                             test = dat$.fid[loo])
    }
    rsmp[[i]] <- as.resampling(rsmp[[i]])
  }
  names(rsmp) <- as.character(repetition)
  as.represampling(rsmp)
}

# parti <- partition_loo_fac_strat(d, repetition = 1:3, ndisc = 5, strat_fac = "croptype", group_fac = "field", n_per_group = 25, seed1 = 123)
# plot(parti, data = d)


# Partitioning method for the Maipo case study:
# - Sampling takes place at the field level
# - A certain number of fields from each crop type is sampled
# - Fields are then partitioned randomly
# This is field-level CV with a small adaptation to
# use only a fraction of fields from the large number of
# fields available in the dataset.
partition_cv_fac_strat <- function(data, 
                                   coords = c("x", "y"), nfold = 10, 
                                   repetition = 1:1,
                                   strat_fac = "class", 
                                   group_fac = "field",
                                   ngroups_per_strat,
                                   relevel = TRUE,
                                   seed1 = NULL) 
{
  if (length(repetition) < 2) {
    repetition <- seq_len(repetition)
  }
  data$.fid <- seq_len(nrow(data))
  rsmp <- list()
  for (i in repetition) {
    if (!is.null(seed1)) {
      set.seed(seed1 + i)
    }
    dat <- resample_fac_strat(data = data, 
                              param = list(strat_fac = strat_fac, group_fac = group_fac, 
                                           ngroups_per_strat = ngroups_per_strat, relevel = relevel))
    resampler <- sample(rep(sample(nfold), length = nrow(dat)), 
                        size = nrow(dat))
    rsmp[[i]] <- list()
    for (j in 1:nfold) {
      sel <- resampler == j
      rsmp[[i]][[j]] <- list(train = dat$.fid[!sel],
                             test  = dat$.fid[sel])
    }
    rsmp[[i]] <- as.resampling(rsmp[[i]])
  }
  names(rsmp) <- as.character(repetition)
  as.represampling(rsmp)
}

# parti <- partition_cv_fac_strat(d, repetition = 1:3, nfold = 4,
#                                 strat_fac = "croptype",
#                                 group_fac = "field",
#                                 ngroups_per_strat = 25,
#                                 seed1 = 123)
# plot(parti, data = d)


# Partitioning method for the Maipo case study:
# - Sampling takes place at the field level
# - A specific number of fields from each crop type is sampled
# - These fields are then spatially partitioned using the k-means
#   method for spatial cross-validation.
partition_kmeans_fac_strat <- function(data, 
                                       coords = c("x", "y"), 
                                       nfold = 10, 
                                       repetition = 1:1,
                                       strat_fac = "class", 
                                       group_fac = "field",
                                       ngroups_per_strat,
                                       relevel = TRUE,
                                       balancing_steps = 1,
                                       order_clusters = TRUE,
                                       seed1 = NULL) 
{
  if (length(repetition) < 2) {
    repetition <- seq_len(repetition)
  }
  balancing_steps <- max(1, balancing_steps)
  
  kmgini <- function(x) {
    p <- x$size/sum(x$size)
    return(1 - sum(p^2))
  }
  
  data$.fid <- seq_len(nrow(data))
  rsmp <- list()
  for (i in repetition) {
    if (!is.null(seed1)) {
      set.seed(seed1 + i)
    }
    dat <- resample_fac_strat(data = data, 
                              param = list(strat_fac = strat_fac, group_fac = group_fac, 
                                           ngroups_per_strat = ngroups_per_strat, relevel = relevel))
    kms <- list()
    for (ik in 1:balancing_steps) {
      kms[[ik]] <- kmeans(dat[, coords], centers = nfold) # additional arguments not supported
    }
    km <- kms[[which.max(sapply(kms, kmgini))]]
    if (order_clusters) {
      o <- rank(km$centers[, 1], ties.method = "first")
      km$cluster <- o[km$cluster]
    }
    resampler <- km$cluster
    rsmp[[i]] <- list()
    for (j in 1:nfold) {
      sel <- resampler == j
      rsmp[[i]][[j]] <- list(train = dat$.fid[!sel],
                             test  = dat$.fid[sel])
    }
    rsmp[[i]] <- as.resampling(rsmp[[i]])
  }
  names(rsmp) <- as.character(repetition)
  as.represampling(rsmp)
}

# parti <- partition_kmeans_fac_strat(d, repetition = 1:3, nfold = 4,
#                                 strat_fac = "croptype",
#                                 group_fac = "field",
#                                 ngroups_per_strat = 25,
#                                 seed1 = 123,
#                                 balancing_steps = 5)
# plot(parti, data = d)


# Partitioning method for the Maipo case study:
# - Resample fields (not grid cells); these are the 'groups'
# - Ensure that within each stratum, exactly ngroups_per_strat
#   fields are sampled. Thus, it can be guaranteed that 
#   e.g. 50 fields of each crop type are included in each 
#   sample.
partition_factor_cv_fac_strat <- function(data, 
                                          coords = c("x", "y"), 
                                          nfold = 10, 
                                          repetition = 1:1,
                                          strat_fac = "class", 
                                          group_fac = "field",
                                          ngroups_per_strat,
                                          relevel = TRUE,
                                          seed1 = NULL) 
{
  if (length(repetition) < 2) {
    repetition <- seq_len(repetition)
  }
  
  data$.fid <- seq_len(nrow(data))
  rsmp <- list()
  for (i in repetition) {
    if (!is.null(seed1)) {
      set.seed(seed1 + i)
    }
    dat <- resample_fac_strat(data = data, 
                              param = list(strat_fac = strat_fac, group_fac = group_fac, 
                                           ngroups_per_strat = ngroups_per_strat, relevel = relevel))
    
    fac <- dat[, group_fac]
    resampler <- sample(rep(sample(nfold), length = nlevels(fac)), 
                        size = nlevels(fac))
    names(resampler) <- levels(fac)
    resampler <- resampler[fac]
    
    rsmp[[i]] <- list()
    for (j in 1:nfold) {
      sel <- resampler == j
      rsmp[[i]][[j]] <- list(train = dat$.fid[!sel],
                             test  = dat$.fid[sel])
    }
    rsmp[[i]] <- as.resampling(rsmp[[i]])
  }
  names(rsmp) <- as.character(repetition)
  as.represampling(rsmp)
}


# Modified summary method with na.rm argument:
summary_rep <- function(object, level = 0, na.rm = TRUE, ...) {
  class(object) <- NULL
  object <- as.data.frame(object)
  if (level <= 0) {
    object <- data.frame(
      mean = sapply(object, mean, na.rm = na.rm),
      sd = sapply(object, sd, na.rm = na.rm),
      median = sapply(object, median, na.rm = na.rm),
      IQR = sapply(object, IQR, na.rm = na.rm)
    )
  }
  object
}




####################################################
### Functions for extracting info from sperrorest objects
####################################################

errdist <- function(x, type = "test") {
  if (is.atomic(x)) {
    return(list(err = numeric(), dist = numeric()))
  } else if (is.null(names(x[[type]]))) {
    return(list(err = numeric(), dist = numeric()))
  } else {
    return(list(err = x[[type]]$bias, dist = x$distance))
  }
}
errdists <- function(x, type = "test") x %>% map(errdist, type = type) %>% dplyr::bind_rows() %>% as.data.frame()


merdist <- function(x) {
  if (is.atomic(x)) {
    return(list(err = numeric(), dist = numeric()))
  } else if (is.null(names(x$test))) {
    return(list(err = numeric(), dist = numeric()))
  } else {
    return(list(err = x$test$error, dist = x$distance))
  }
}
merdists <- function(x) x %>% map(merdist) %>% dplyr::bind_rows() %>% as.data.frame()


biasdist <- function(x) {
  if (is.atomic(x)) {
    return(list(err = numeric(), dist = numeric()))
  } else if (is.null(names(x$test))) {
    return(list(err = numeric(), dist = numeric()))
  } else {
    return(list(err = x$test$bias, dist = x$distance))
  }
}
biasdists <- function(x) x %>% map(biasdist) %>% dplyr::bind_rows() %>% as.data.frame()



impdist <- function(x) {
  if (is.atomic(x)) {
    return(NULL)
  } else {
    imp <- x$error
    names(imp) <- rownames(x)
    return(as.list(imp))
  }
}
impdists <- function(x) x %>% map(impdist) %>% dplyr::bind_rows() %>% as.data.frame()

impdist2 <- function(x) {
  if (is.atomic(x)) {
    return(NULL)
  } else {
    imp <- x$bias
    names(imp) <- rownames(x)
    return(as.list(imp))
  }
}
impdists2 <- function(x) x %>% map(impdist2) %>% dplyr::bind_rows() %>% as.data.frame()


mers <- function(x, which = "err", breaks = 25,
                 smoother = "MA", # or loess
                 span = 0.5, trafo = "sqrt", logy = FALSE,
                 negate = FALSE) {
  if (length(breaks) <= 1)
    breaks <- seq(min(x$dist, na.rm = TRUE), 
                  max(x$dist, na.rm = TRUE), 
                  length = breaks)
  if (smoother == "loess") {
    if (!is.null(trafo)) {
      if (trafo == "log") {
        x$dist <- log10(x$dist)
        if (length(breaks) > 1) breaks <- log10(breaks)
      } else if (trafo == "sqrt") {
        x$dist <- sqrt(x$dist)
        if (length(breaks) > 1) breaks <- sqrt(breaks)
      }
    }
    dat <- data.frame(err = x[,which], dist = x$dist)  
    if (negate) dat$err <- dat$err * (-1)
    if (logy) dat$err <- log10(dat$err+0.001)
    lo <- loess(err ~ dist, data = dat, span = span)
    mer <- predict(lo, newdata = data.frame(dist = breaks))
    if (logy) mer <- exp(mer) - 0.001
    res <- data.frame(
      dist = breaks,
      mer = mer
    )
    if (!is.null(trafo)) {
      if (trafo == "log") res$dist <- 10^res$dist
      if (trafo == "sqrt") res$dist <- res$dist^2
    }
  } else {
    intv <- cut(x$dist, breaks = c(0,breaks))
    mer <- unname(tapply(x[,which], intv, mean))
    breaks <- unname(tapply(x$dist, intv, median))
    res <- data.frame(
      dist = breaks,
      mer = mer
    )
    if (negate) res$mer <- res$mer * (-1)
    res$mer <- smth(res)$smth_mer
  }
  res
}


# Function for calculating RMSE within bins:
# The x column identified by the 'which' argument
# contains the e_i := y_i - y^_i error values.
rmses0 <- function(x, which = "err", breaks = 25) {
  intv <- cut(x$dist, breaks = breaks)
  rmse <- function(lev) sqrt(mean((x[intv==lev, which])^2, na.rm = TRUE))
  rmses <- levels(intv) %>% map(rmse) %>% unlist()
  data.frame(
    dist = unname(tapply(x$dist, intv, median)),
    rmse = unname(rmses)
  )
}

# Function for calculating RMSE within bins;
# here, the x column identified by the 'which' argument
# contains the predictions, and the 'obs' column the
# observations.
rmses <- function(x, which = "err", breaks = 25) {
  intv <- cut(x$dist, breaks = breaks)
  rmse <- function(lev) sqrt(mean((x[intv==lev, which] - x[intv==lev,"obs"])^2, na.rm = TRUE))
  rmses <- levels(intv) %>% map(rmse)
  data.frame(
    dist = unname(tapply(x$dist, intv, median)),
    rmse = unlist(rmses)
  )
}


# Simple smoothing function based on a weighted moving average:
smth <- function(x, exclude = "dist") {
  nms <- colnames(x)
  nms <- nms[ !(nms %in% exclude) ]
  for (nm in nms) {
    x$TEMPNAME <- 0.5 * x[,nm] + 
      0.25 * c(x[-1, nm], x[nrow(x), nm]) + 
      0.25 * c(x[1, nm], x[-nrow(x), nm])
    colnames(x)[ncol(x)] <- paste0("smth_", nm)
  }
  x
}

####################################################
### Error functions
####################################################


# Simplified error function for Maipo case study:
err_maipo <- function(obs, pred) {
  acc <- mean(as.double(obs == pred))
  list(
    accuracy = acc,
    error = 1.0 - acc
  )
}

# Simplified error function for Meuse case study:
err_meuse <- function(obs, pred) {
  list(
    bias = mean(obs - pred),
    obs = mean(obs),
    pred = mean(pred),
    count = length(obs)
  )
}



####################################################
### Models
####################################################


# Ordinary kriging, implemented as a wrapper that ignores
# any external drift variables in the model formula:
mykrigeOK <- function(formula, data, locations = ~x + y, ...) {
  # Create formula of type y ~ 1:
  OKformula <- as.formula(paste(all.vars(formula)[1], "~ 1"))
  alexmisc::mykrige(formula = OKformula, data = data, locations = locations, ...)
}

# mygwr model, i.e. a wrapper for GWR that automatically
# tunes its bandwidth parameter.
# x/y coordinates are only used as spatial reference, not
# as predictors themselves since this would make the model
# too flexible and probably nearly unidentifiable.
mygwr <- function(formula, data, coords = ~x+y, crs = CRS("+init=epsg:28992"),
                  remove_coords_from_formula = TRUE,
                  verbose = 0) {
  if (remove_coords_from_formula) {
    xvars <- all.vars(formula)[-1]
    coordvars <- all.vars(coords)
    xvars <- xvars[ !(xvars %in% coordvars) ]
    formula <- as.formula(paste(all.vars(formula)[1], "~",
                                paste(xvars, collapse = "+")))
  }  
  if (verbose >= 1)
    cat("GWR: n = ", nrow(data), ", bandwidth = ", sep = "")
  coordinates(data) <- coords
  proj4string(data) <- crs
  fitbw <- spgwr::gwr.sel(formula = formula, data = data, 
                          adapt = FALSE, verbose = verbose > 1)
  if (verbose >= 1)
    cat(fitbw, ", ", sep = "")
  fitgwr <- spgwr::gwr(formula = formula, data = data, bandwidth = fitbw)
  if (verbose >= 1)
    cat("model has been fitted.\n", sep = "")
  res <- list(
    bandwidth = fitbw,
    fit = fitgwr,
    formula = formula,
    data = data,
    coords = coords,
    crs = crs
  )
  class(res) <- "mygwr"
  res
}

# Predict method for mygwr objects:
predict.mygwr <- function(object, newdata, verbose = 0) {
  if (verbose >= 1)
    cat("Predict from GWR: n_train = ", nrow(object$data), 
        ", n_pred = ", nrow(newdata), ", ", sep = "")
  pred <- c()
  if (nrow(newdata) > 0) {
    coordinates(newdata) <- object$coords
    proj4string(newdata) <- object$crs
    predgwr <- spgwr::gwr(formula = object$formula, 
                          data = object$data, predictions = TRUE, 
                          bandwidth = object$bandwidth, 
                          fit.points = newdata)
    pred <- predgwr$SDF$pred
  }
  if (verbose >= 1)
    cat("predicted", length(pred), "values.\n")
  pred
}



# 1-nearest-neighbour method:
myNN <- function(data, newdata, response, coords = c("x", "y"), 
                 log_dist = FALSE, min_dist = 0.1) {
  nn1 <- function(i) {
    di <- sqrt( (data[, coords[1]] - newdata[i, coords[1]])^2 +
                (data[, coords[2]] - newdata[i, coords[2]])^2 )
    wh <- nnet::which.is.max(-di)
    if (length(wh) != 1) {
      cat("data:\n\n")
      print(data)
      cat("\n\n\nnewdata:\n\n")
      print(newdata)
    }
    stopifnot(length(wh) == 1)
    data.frame(nn_pred = data[wh, response], nn_dist = di[wh])
  }
  if (nrow(newdata) == 0) 
    return(data.frame(nn_pred = numeric(), nn_dist = numeric()))
  
  res <- dplyr::bind_rows(lapply(1:nrow(newdata), nn1))
  stopifnot(all(!is.na(res$nn_pred)))
  
  if (log_dist) {
    res$nn_dist[ res$nn_dist == 0 ] <- min_dist
    res$nn_dist <- log10(res$nn_dist)
  }
  res
}

# 1-nearest neighbour classification:
# Fitting function simply returns the data and metadata.
nn <- function(formula, data, coords = c("utmx","utmy"), ...) {
  response <- all.vars(formula)[1]
  fit <- list(data = na.omit(data[, c(coords, response)]), 
              response = response,
              coords = coords) 
  class(fit) <- "nn"
  fit
}

# Predict method for my implementation of 1-nearest-neighbour classification:
predict.nn <- function(object, newdata) {
  myNN(data = object$data, newdata = newdata, 
       response = object$response, coords = object$coords)$nn_pred
}

# splda implements the hybrid NN-LDA method:
splda <- function(formula, data, coords = c("x", "y"),
                  maxdist = 100, remove_xy = TRUE, ...)
{
  yvar <- all.vars(formula)[1]
  if (remove_xy) {
    vnms <- all.vars(formula)[-1]
    vnms <- vnms[ !(vnms %in% coords) ]
    formula <- as.formula(paste(yvar, "~",
                                paste(vnms, collapse = "+")))
  }
  fit <- MASS::lda(formula = formula, data = data, ...)
  fit <- list(
    ldafit = fit,
    formula = formula,
    xy = data[ , coords],
    response = data[, yvar],
    coords = coords,
    maxdist = maxdist
  )
  class(fit) <- "splda"
  fit
}

# Predict method for fitted splda objects, i.e. NN-LDA models:
predict.splda <- function(object, newdata) {
  mynn <- function(i) {
    di <- sqrt((object$xy[,1] - newxy[i, 1])^2 +
               (object$xy[,2] - newxy[i, 2])^2)
    nnpred <- object$response[nnet::which.is.max(-di)]
    list(dist = min(di),
         nn = nnpred)
  }
  newxy <- newdata[, object$coords]
  nnpred <- seq_len(nrow(newdata)) %>% map(mynn) %>% bind_rows() %>% as.data.frame()
  pred <- nnpred$nn
  sel <- nnpred$dist <= object$maxdist
  if (!all(sel))
    pred[!sel] <- predict(object = object$ldafit, newdata = newdata[!sel,])$class
  pred
}

# Predict method for fitted warped splda object;
# this is probably only needed because predict.splda
# as defined here in this R script (and not in a package)
# is not an officially registered S3 method for splda.
predict.wsplda <- function(object, newdata) {
  newdata <- unwarp(newdata, warper = object$warper)
  predict.splda(object$fit, newdata = newdata)
}



# Hybrid ordinary kriging - random forest model
okrf <- function(formula, data, 
                  taper_radius = 500,
                  locations = ~x+y, 
                  svgm = "Sph", range = NA,
                  nsratio = 0.25, cressie = TRUE,
                  fixed = FALSE, fit.ranges = TRUE, nmax = Inf,
                  ...) {
  # Fit ordinary kriging model:
  okfit <- mykrigeOK(formula = formula, data = data, locations = locations,
                     svgm = svgm, range = range, nsratio = nsratio,
                     cressie = cressie, fixed = fixed, fit.ranges = fit.ranges,
                     nmax = nmax)
  # Fit random forest model:
  rffit <- randomForest::randomForest(formula = formula, data = data, ...)
  # Return fitted models along with data and other arguments
  # that are needed to taper both models in the prediction step:
  obj <- list(
    okfit = okfit, 
    rffit = rffit,
    data = data,
    coords = all.vars(locations),
    taper_radius = taper_radius
  )
  class(obj) <- "okrf"
  obj
}

# Predict method for fitted okrf models (OK-RF):
predict.okrf <- function(object, newdata, ...) {
  mindi <- rep(NA, nrow(newdata))
  # Calculate each prediction location's minimum distance
  # to training sample:
  for (i in seq_len(nrow(newdata))) {
    mindi[i] <- min(sqrt(
      (object$data[, object$coords[1]] - newdata[i, object$coords[1]])^2 + 
        (object$data[, object$coords[2]] - newdata[i, object$coords[2]])^2
    ))
  }
  okpred <- predict(object = object$okfit, newdata = newdata)
  rfpred <- predict(object = object$rffit, newdata = newdata)
  wgt <- mindi / object$taper_radius
  wgt[wgt > 1] <- 1
  pred <- wgt * rfpred + (1-wgt) * okpred
  pred
}


# Predict function for LDA that predicts class membership:
predict_lda <- function(object, newdata)
  predict(object, newdata)$class
