####################################################
### Spatial model diagnostics: Meuse
####################################################
### Plot results from Meuse case study
####################################################
### (c) 2021 Alexander Brenning
####################################################

# Load packages:
library("sperrorest")
library("raster")

# Character expansion factor for text:
CEX <- 0.75

# meuse results:
MEUSE <- readRDS("meuse_smry.rds")


####################################################
### Figure 1
####################################################
### Spatial prediction maps of $logZn$ in the Meuse 
### case study using four selected geostatistical, 
### ML, and hybrid models.
####################################################

par(mfrow = c(2,2), cex = CEX, cex.main = 1,
    mar = c(3.5, 3, 3, 3))
rg <- c(min(c(MEUSE$spnewdata$pred_ked, MEUSE$spnewdata$pred_ok, MEUSE$spnewdata$pred_rf, MEUSE$spnewdata$pred_okrf)),
        max(c(MEUSE$spnewdata$pred_ked, MEUSE$spnewdata$pred_ok, MEUSE$spnewdata$pred_rf, MEUSE$spnewdata$pred_okrf)))
brks <- pretty(rg, n = 10)
cols <- hcl.colors(length(brks) - 1)
plot(MEUSE$r_pred_ok, main = "Ordinary kriging", breaks = brks, col = cols)
plot(MEUSE$r_pred_ked, main = "Kriging with external drift", breaks = brks, col = cols)
plot(MEUSE$r_pred_okrf, main = "Hybrid OK-RF", breaks = brks, col = cols)
plot(MEUSE$r_pred_rf, main = "Random forest", breaks = brks, col = cols)


####################################################
### Figure 2
####################################################
### Spatial prediction error profiles for the 
### prediction of logZn in the Meuse case study. 
### The x axes are square-root transformed. 
### LOO-CV and other point estimators are placed 
### close to their mean prediction distance. 
### Left: All models; right: detailed view of KED, 
### MLR, RF, and OK--RF. Models: NN (brown), 
### OK (grey), KED (black), MLR (blue), GWR (light 
### blue), RF (dark green), OK--RF (light green).
####################################################

adj0 <- 0.125
xadj <- list(
  NN = 0,
  LM = 0,
  RF = 0,
  OK = -adj0,
  OKRF = adj0,
  KED = adj0,
  GWR = 0
)
xlim <- c(20,1000)
ylim <- c(0.0, 0.57)
xs <- c(0, 60, 200, 400, 700, 1000)
xlabels <- as.character(xs)
par(mfrow = c(1,2), mar = c(3.8,3.8,0.5,0.7), mgp = c(2,.7,0), cex = CEX)
plot(NA, xlim = sqrt(xlim), ylim = ylim, 
     xlab = "Prediction distance [m]", ylab = "RMSE", 
     xaxt = "n")
axis(1, at = sqrt(xs), label = xlabels)
for (p in MEUSE$pal) {
  lines(smth_rmse ~ sqrt(dist), MEUSE$rmse[[p$name]], col = p$col, lwd = 2, lty = p$lty)
  points(sqrt(22)+xadj[[p$name]], MEUSE$rmse_train[[p$name]], col = p$col, pch = 19)
  points(sqrt(MEUSE$mean_mindist-10)+xadj[[p$name]], MEUSE$rmse_loo[[p$name]], col = p$col, pch = 19)
  points(sqrt(MEUSE$mean_mindist+10)+xadj[[p$name]], MEUSE$rmse_cv[[p$name]], col = p$col, pch = 19)
  points(sqrt(MEUSE$mean_spcvdist)+xadj[[p$name]], MEUSE$rmse_spcv[[p$name]], col = p$col, pch = 19)
}
ytxt <- ylim[2] - 0.01
text(sqrt(22), ytxt, "Resubstitution", col = "black", pos = 4, srt = 270, offset = 0, adj = 0.5)
text(sqrt(MEUSE$mean_mindist-15), ytxt, "LOO-CV", col = "black", pos = 4, srt = 270, offset = 0, adj = 0.5)
text(sqrt(MEUSE$mean_cvdist+15), ytxt, "Random CV", col = "black", pos = 4, srt = 270, offset = 0, adj = 0.5)
text(sqrt(MEUSE$mean_spcvdist), ytxt, "Spatial\nCV", col = "black", pos = 4, srt = 270, offset = 0, adj = 0.5)
text(sqrt(650), 0.43, "NN", col = MEUSE$pal$NN$col)
text(sqrt(650), 0.35, "OK", col = MEUSE$pal$OK$col)
#text(sqrt(670), 0.27, "GWR", pos = 4, col = MEUSE$pal$GWR$col)



ylim <- c(0.05, 0.25)
plot(NA, xlim = sqrt(xlim), ylim = ylim, xlab = "Prediction distance [m]", ylab = "RMSE", xaxt = "n")
axis(1, at = sqrt(xs)[-1], label = xlabels[-1])
for (p in MEUSE$pal[c("LM", "KED", "OKRF", "RF")]) {
  lines(smth_rmse ~ sqrt(dist), MEUSE$rmse[[p$name]], col = p$col, lwd = 2, lty = p$lty)
  points(sqrt(22)+xadj[[p$name]], MEUSE$rmse_train[[p$name]], col = p$col, pch = 19)
  points(sqrt(MEUSE$mean_mindist-10)+xadj[[p$name]], MEUSE$rmse_loo[[p$name]], col = p$col, pch = 19)
  points(sqrt(MEUSE$mean_mindist+10)+xadj[[p$name]], MEUSE$rmse_cv[[p$name]], col = p$col, pch = 19)
  points(sqrt(MEUSE$mean_spcvdist)+xadj[[p$name]], MEUSE$rmse_spcv[[p$name]], col = p$col, pch = 19)
}
ytxt <- ylim[2]
text(sqrt(22), ytxt, "Resubstitution", col = "black", pos = 4, srt = 270, offset = 0, adj = 0.5)
text(sqrt(MEUSE$mean_mindist-15), ytxt, "LOO-CV", col = "black", pos = 4, srt = 270, offset = 0, adj = 0.5)
text(sqrt(MEUSE$mean_mindist+15), ytxt, "Rnd. CV", col = "black", pos = 4, srt = 270, offset = 0, adj = 0.5)
text(sqrt(MEUSE$mean_spcvdist), ytxt, "Spatial\nCV", col = "black", pos = 4, srt = 270, offset = 0, adj = 0.5)
text(sqrt(60), 0.15, "MLR", col = MEUSE$pal$LM$col)
text(sqrt(60), 0.075, "KED", col = MEUSE$pal$KED$col, pos = 4)
text(sqrt(110), 0.12, "RF", pos = 4, col = MEUSE$pal$RF$col)
text(sqrt(40), 0.127, "OK-RF", col = MEUSE$pal$OKRF$col)



####################################################
### Figure 3
####################################################
### Spatial variable importance profiles in the 
### prediction of $logZn$ in the Meuse case study. 
### Models: NN (brown), OK (grey), KED (black), 
### MLR (blue), GWR (light blue), RF (dark green), 
### OK--RF (light green).
####################################################


xlim <- c(60,1000)
par(mfrow = c(2,2), mar = c(3.8,3.8,2.5,0.7), mgp = c(2,.7,0), cex = CEX)

vnms <- stringr::str_subset(colnames(MEUSE$dmer$LM), "^smth_")
for (vnm in vnms) {
  dylim <- c(0.0, 0.45)
  if (any(vnm == c("smth_sqrt.dist", "smth_elev")))
    dylim <- c(0, 0.22)
  
  the_vnm <- gsub("smth_", "", vnm)
  
  plot(NA, xlim = sqrt(xlim), ylim = dylim, xaxt = "n",
       xlab = "Prediction distance [m]", ylab = "Increase in RMSE",
       main = the_vnm)
  xs <- c(70, 200, 400, 700, 1000)
  axis(side = 1, at = sqrt(xs), labels = xs)
  for (p in MEUSE$pal) {
    lines(sqrt(MEUSE$dmer[[p$name]]$dist), 
          MEUSE$dmer[[p$name]][, vnm],
          lty = p$lty, col = p$col, lwd = 2)
    points(sqrt(MEUSE$mean_mindist-10)+xadj[[p$name]], -MEUSE$imp_loo[p$name, the_vnm], col = p$col, pch = 19)
    points(sqrt(MEUSE$mean_mindist+10)+xadj[[p$name]], -MEUSE$imp_cv[p$name, the_vnm], col = p$col, pch = 19)
    points(sqrt(MEUSE$mean_spcvdist)+xadj[[p$name]], -MEUSE$imp_spcv[p$name, the_vnm], col = p$col, pch = 19)
  }
  if (the_vnm == "elev") {
    ytxt <- dylim[2]
    text(sqrt(MEUSE$mean_mindist-15), ytxt, "LOO-CV", col = "black", pos = 4, srt = 270, offset = 0, adj = 0.5)
    text(sqrt(MEUSE$mean_mindist+15), ytxt, "Random CV", col = "black", pos = 4, srt = 270, offset = 0, adj = 0.5)
    text(sqrt(MEUSE$mean_spcvdist), ytxt, "Spatial CV", col = "black", pos = 4, srt = 270, offset = 0, adj = 0.5)
  }
}
