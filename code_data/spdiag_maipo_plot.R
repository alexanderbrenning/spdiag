####################################################
### Spatial model diagnostics: Maipo
####################################################
### Plot results from Maipo case study
####################################################
### (c) 2021 Alexander Brenning
####################################################

# Load packages:
library("sperrorest")
library("raster")

# Character expansion factor for text:
CEX <- 0.75

# maipo results:
MAIPO <- readRDS("maipo_smry.rds")

library("stringr")
MAIPO$nbands <- length(stringr::str_subset(MAIPO$xvars, "^b"))
MAIPO$nindices <- length(stringr::str_subset(MAIPO$xvars, "^nd"))


####################################################
### Figure 4
####################################################
### Spatial prediction error profiles for the
### classification of crop type in the Maipo case 
### study. Point estimates from different CV types 
### are plotted at their mean prediction distances. 
### Results of LOO-CV and non-spatial random CV are 
### visually indistinguishable.
####################################################

xlim <- c(0, 10000) #MAXDIST)
xs <- c(0, 250, 1000, 2500, 5000, 10000) #, 15000)
ylim <- c(0,0.22)
par(mfrow = c(1,1), mar = c(3.8,3.8,0.5,0.7), mgp = c(2,.7,0), cex = CEX)
plot(mer ~ sqrt(dist), MAIPO$mer_splda[MAIPO$mer_splda$dist < 500,], type = "l", 
     col = MAIPO$pal$SPLDA$col, lty = MAIPO$pal$SPLDA$lty, lwd = 2, 
     xlim = sqrt(xlim), ylim = ylim, 
     xlab = "Prediction distance [m]", ylab = "Error rate", 
     xaxt = "n")
axis(1, at = sqrt(xs), label = xs)
lines(mer ~ sqrt(dist), MAIPO$mer_rf[MAIPO$mer_rf$dist < MAIPO$maxdist+200,,], 
      col = MAIPO$pal$RF$col, lwd = 2, lty = MAIPO$pal$RF$lty)
lines(mer ~ sqrt(dist), MAIPO$mer_lda[MAIPO$mer_lda$dist < MAIPO$maxdist+200,,], 
      col = MAIPO$pal$LDA$col, lwd = 2, lty = MAIPO$pal$LDA$lty)

points(sqrt(rep(MAIPO$mean_spcvdist, 2)), c(MAIPO$mer_lda_spcv, MAIPO$mer_rf_spcv), 
       pch = 19, col = c(MAIPO$pal$LDA$col, MAIPO$pal$RF$col))
points(sqrt(rep(MAIPO$mean_cvdist, 2)), c(MAIPO$mer_lda_cv, MAIPO$mer_rf_cv), 
       pch = 19, col = c(MAIPO$pal$LDA$col, MAIPO$pal$RF$col))
# Visually indistinguishable from LOO-CV:
# points(sqrt(rep(MAIPO$mean_ncvdist, 3)), c(MAIPO$mer_lda_ncv, MAIPO$mer_rf_ncv, MAIPO$mer_sploo_ncv), 
#        pch = 19, col = c(MAIPO$pal$LDA$col, MAIPO$pal$RF$col, MAIPO$pal$SPLDA$col))
points(sqrt(c(30, 30, 30)), c(MAIPO$mer_lda_loo, MAIPO$mer_rf_loo, MAIPO$mer_splda_loo), 
       pch = 19, col = c(MAIPO$pal$LDA$col, MAIPO$pal$RF$col, MAIPO$pal$SPLDA$col))
points(c(0,-0.25,0.25), 
       c(MAIPO$mer_lda_train, MAIPO$mer_rf_train, MAIPO$mer_splda_train), 
       pch = 16, col = c(MAIPO$pal$LDA$col, MAIPO$pal$RF$col, MAIPO$pal$SPLDA$col))

text(sqrt(2500), 0.09, "LDA", col = MAIPO$pal$LDA$col, pos = 1)
text(sqrt(2500), 0.165, "Random forest", col = MAIPO$pal$RF$col, pos = 3)
text(sqrt(100), 0.02, "NN-LDA", col = MAIPO$pal$SPLDA$col, pos = 4)

text(sqrt(0), 0.22, "Resubstitution", col = "black", pos = 4, srt = 270, offset = 0, adj = 0.5)
# text(sqrt(MAIPO$mean_ncvdist), 0.22, "Random CV", col = "black", pos = 4, srt = 270, offset = 0, adj = 0.5)
# text(sqrt(30), 0.22, "LOO-CV", col = "black", pos = 4, srt = 270, offset = 0, adj = 0.5)
text(sqrt(30), 0.22, "LOO-CV, random CV", col = "black", pos = 4, srt = 270, offset = 0, adj = 0.5)
text(sqrt(MAIPO$mean_cvdist), 0.22, "Field-level\nCV", col = "black", pos = 1)
text(sqrt(MAIPO$mean_spcvdist), 0.22, "Spatial\nCV", col = "black", pos = 1)



####################################################
### Figure 5
####################################################
### Spatial variable importance profiles in crop
### classification in the Maipo case study for the 
### first PC of early- and mid-season features and 
### the first two PCs of late-season features. 
### Models: RF (green), LDA (blue), NN--LDA 
### (light blue).
####################################################

dylim <- c(0, 0.4) 

MAIPO$imp_names_plot <- c("Early1", "Mid1", "Late1", "Late2")
nrw <- ifelse(length(MAIPO$imp_names_plot) <= 6, 2, 3)
ncl <- ceiling(length(MAIPO$imp_names_plot) / nrw)
par(mfrow = c(nrw, ncl), mar = c(3.8,3.8,2.5,0.7), mgp = c(2,.7,0), cex = CEX)

for (vnm in MAIPO$imp_names_plot) {
  plot(NA, xlim = sqrt(xlim), ylim = dylim,
       xlab = "Prediction distance [m]", ylab = "Decrease in accuracy",
       main = vnm,
       xaxt = "n")
  axis(1, at = sqrt(xs), label = xs)
  
  for (m in names(MAIPO$dmer)) {
    lines(sqrt(MAIPO$dmer[[m]]$dist), MAIPO$dmer[[m]][, vnm],
          col = MAIPO$pal[[m]]$col, 
          lty = MAIPO$pal[[m]]$lty,
          lwd = 2)
  }
}



####################################################
### Appendix
####################################################
### Histograms of prediction distances in the Maipo 
### case study. Left: Distribution based on randomly 
### selecting 100 fields for training and the 
### remaining 300 fields as prediction locations. 
### Center and right: Distributions for field-level 
### CV and $k$-means-based spatial CV, respectively. 
### All histograms are aggregated over 50 samples or 
### CV repetitions.
####################################################

par(mfrow = c(1, 3), mar = c(3.6,3.6,2.7,0.5), mgp = c(2,.7,0), cex = CEX, cex.main = 1.0)
xlim <- c(0, max(MAIPO$spcvdist))
brks <- seq(xlim[1], xlim[2], len = 40)
hist(MAIPO$preddist, breaks = brks, 
     xlim = xlim, col = "lightblue", 
     xlab = "Prediction distance [m]",
     freq = FALSE,
     main = "Prediction in entire area")
text(13000, 4.5e-04, paste0("Mean: ", round(MAIPO$mean_preddist), " m"))
lines(rep(MAIPO$mean_preddist, 2), c(0, 8e-04), lty = "dashed")
hist(MAIPO$cvdist, breaks = brks, 
     xlim = xlim, col = "lightblue", 
     xlab = "Prediction distance [m]",
     freq = FALSE,
     main = "Prediction in field-level CV")
text(13000, 4.3e-04, paste0("Mean: ", round(MAIPO$mean_cvdist), " m"))
lines(rep(MAIPO$mean_cvdist, 2), c(0, 8e-04), lty = "dashed")
hist(MAIPO$spcvdist, breaks = brks, 
     xlim = xlim, col = "lightblue", 
     xlab = "Prediction distance [m]",
     freq = FALSE,
     main = "Prediction in spatial CV")
text(16000, 0.00014, paste0("Mean:\n", round(MAIPO$mean_spcvdist), " m"))
lines(rep(MAIPO$mean_spcvdist, 2), c(0,0.0002), lty = "dashed")

