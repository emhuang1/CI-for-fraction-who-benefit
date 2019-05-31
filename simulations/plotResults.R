rm(list=ls())

##Vector of candidates psi in grid G[0,1] (with point at every hundredth)
psi <- seq(from = 0, to = 1, by = 0.01)


##############################################################################################################
##Functions to generate coverage probability plots##
##############################################################################################################

plotCov_singleSetting <- function(n, lb, ub, nsim, homWd, plotTitle, xregion, yregion){
  setwd(homWd)
  setwd("ourMethod/results")
  load(paste("res",n,".Rdata", sep = ""))
  if (plotTitle == "Setting D"){     
    CI.our <- ci ##CI has previously been computed
  } else {
    CI.our <- matrix(NA, nrow = nsim, ncol = 2) ##compute CI from the confidence set
    for (sim in 1:nsim){
      temp <- which(ci[sim,]==1)
      if (length(temp) >= 1){
        CI.our[sim,] <- c(psi[min(temp)], psi[max(temp)])
      } else{
        CI.our[sim,] <- c(NA,NA)
      }  
    }
  }
  
  setwd(homWd)
  setwd("mOfn/results")
  load(paste("res",n,".Rdata", sep = ""))
  CI.mOfn <- ci[1:nsim,]
  rm(ci)
  
  setwd(homWd)
  setwd("manski_horowitz/results")
  load(paste("res",n,".Rdata", sep = ""))
  CI.horowitz <- ci[1:nsim,]
  rm(ci)
  
  cov.our <- rep(NA,length(psi))
  cov.mOfn.25 <- rep(NA,length(psi))
  cov.mOfn.50 <- rep(NA,length(psi))
  cov.mOfn.75 <- rep(NA,length(psi))
  cov.mOfn.90 <- rep(NA,length(psi))
  cov.mOfn.100 <- rep(NA,length(psi))
  cov.horowitz <- rep(NA, length(psi))
  
  for (i in 1:length(psi)){
    cov.our[i] <- sum(psi[i] >= CI.our[,1] & psi[i] <= CI.our[,2])/nsim
    cov.mOfn.100[i] <- sum(psi[i] >= CI.mOfn[,1] & psi[i] <= CI.mOfn[,2])/nsim
    cov.mOfn.90[i] <- sum(psi[i] >= CI.mOfn[,3] & psi[i] <= CI.mOfn[,4])/nsim
    cov.mOfn.75[i] <- sum(psi[i] >= CI.mOfn[,5] & psi[i] <= CI.mOfn[,6])/nsim
    cov.mOfn.50[i] <- sum(psi[i] >= CI.mOfn[,7] & psi[i] <= CI.mOfn[,8])/nsim
    cov.mOfn.25[i] <- sum(psi[i] >= CI.mOfn[,9] & psi[i] <= CI.mOfn[,10])/nsim
    cov.horowitz[i] <- sum(psi[i] >= CI.horowitz[,1] & psi[i] <= CI.horowitz[,2])/nsim
  }
  
  plot(0,0, type = "n", xlim = xregion, ylim = yregion, xaxt = "n", yaxt = "n", xlab = NA, ylab = NA, yaxs = "i", bty = "l")
  rect(lb,0,ub,1, col = "grey90", border = NA)
  par(new = TRUE)
  plot(psi,cov.mOfn.25, xlim = xregion, ylim = yregion, main = NA, type = "l",lty = 4, xlab = expression(psi), ylab = NA, lwd = 2, cex.lab=1.5, cex.axis=1.5, cex.main=1,yaxs = "i", bty = "l")
  title(ylab = expression(paste("Probability that the CI contains ",psi)), line=2.5, cex.lab=1.5)
  par(new = TRUE)
  ##plot(psi,cov.mOfn.50, xlim = xregion, ylim = yregion, xaxt = "n", yaxt = "n", main = NA, xlab = NA, ylab = NA, type = "l", lty = 5, lwd = 2, cex.lab=1.5, cex.main=1, yaxs = "i", bty = "l")
  ##par(new = TRUE)
  ##plot(psi,cov.mOfn.75, xlim = c(0,1), ylim = c(0,1), xaxt = "n", yaxt = "n", main = NA, xlab = NA, ylab = NA, col = "purple", type = "l", lwd = 1.5, cex.lab=2, cex.main=2, yaxs = "i")
  ##par(new = TRUE)
  ##plot(psi,cov.mOfn.90, xlim = c(0,1), ylim = c(0,1), xaxt = "n", yaxt = "n", main = NA, xlab = NA, ylab = NA, col = "orange", type = "l", lwd = 1.5, cex.lab=2, cex.main=2, yaxs = "i")
  ##par(new = TRUE)
  plot(psi,cov.mOfn.100, xlim = xregion, ylim = yregion, xaxt = "n", yaxt = "n", main = NA, xlab = NA, ylab = NA, type = "l", lty = 3, lwd = 2, cex.lab=1.5, cex.main=1, yaxs = "i", bty = "l")
  par(new = TRUE)
  plot(psi,cov.our, xlim = xregion, ylim = yregion, xaxt = "n", yaxt = "n", main = NA, xlab = NA, ylab = NA, lty = 1, type = "l", lwd = 1.5, cex.lab=1.5, cex.main=1, yaxs = "i", bty = "l")
  par(new = TRUE)
  plot(psi,cov.horowitz, xlim = xregion, ylim = yregion, xaxt = "n", yaxt = "n", main = NA, xlab = NA, ylab = NA, type = "l", lty = 5, lwd = 1.5, cex.lab=1.5, cex.main=1, yaxs = "i", bty = "l")
  abline(h = 0.95, lty = 2)
  
  # if (plotTitle == "Setting D"){
  #  legend("bottomleft", c("our method","m = n", "m = 0.25n","Horowitz-Manski"), lty = c(1,3,4,5), lwd = c(1.5,2,2,1.5))
  # } else {
  #  legend("bottomright", c("our method","m = n", "m = 0.25n","Horowitz-Manski"), lty = c(1,3,4,5), lwd = c(1.5,2,2,1.5))
  # }
  ##UNCOMMENT LINES ABOVE FOR MAIN PAPER PLOTS
  
  if (plotTitle == "Setting D"){
    legend("left", c("our method","m = n", "m = 0.25n","Horowitz-Manski"), lty = c(1,3,4,5), lwd = c(1.5,2,2,1.5), cex = 0.8)
  } else {
    legend("right", c("our method","m = n", "m = 0.25n","Horowitz-Manski"), lty = c(1,3,4,5), lwd = c(1.5,2,2,1.5), cex = 0.8)
  }
  ##UNCOMMENT LINES ABOVE FOR SUPP MAT PLOTS
}

plotCov_singleSetting_zoomIn <- function(n, lb, ub, nsim, homWd, plotTitle, xregion, yregion){
  setwd(homWd)
  setwd("ourMethod/results")
  load(paste("res",n,".Rdata", sep = ""))
  if (plotTitle == "Setting D"){
    CI.our <- ci ##CI has previously been computed
  } else {
    CI.our <- matrix(NA, nrow = nsim, ncol = 2) ##compute CI from the confidence set
    for (sim in 1:nsim){
      temp <- which(ci[sim,]==1)
      if (length(temp) >= 1){
        CI.our[sim,] <- c(psi[min(temp)], psi[max(temp)])
      } else{
        CI.our[sim,] <- c(NA,NA)
      }  
    }
  }
  
  setwd(homWd)
  setwd("mOfn/results")
  load(paste("res",n,".Rdata", sep = ""))
  CI.mOfn <- ci[1:nsim,]
  rm(ci)
  
  setwd(homWd)
  setwd("manski_horowitz/results")
  load(paste("res",n,".Rdata", sep = ""))
  CI.horowitz <- ci[1:nsim,]
  rm(ci)
  
  cov.our <- rep(NA,length(psi))
  cov.mOfn.25 <- rep(NA,length(psi))
  cov.mOfn.50 <- rep(NA,length(psi))
  cov.mOfn.75 <- rep(NA,length(psi))
  cov.mOfn.90 <- rep(NA,length(psi))
  cov.mOfn.100 <- rep(NA,length(psi))
  cov.horowitz <- rep(NA, length(psi))
  for (i in 1:length(psi)){
    cov.our[i] <- sum(psi[i] >= CI.our[,1] & psi[i] <= CI.our[,2])/nsim
    cov.mOfn.100[i] <- sum(psi[i] >= CI.mOfn[,1] & psi[i] <= CI.mOfn[,2])/nsim
    cov.mOfn.90[i] <- sum(psi[i] >= CI.mOfn[,3] & psi[i] <= CI.mOfn[,4])/nsim
    cov.mOfn.75[i] <- sum(psi[i] >= CI.mOfn[,5] & psi[i] <= CI.mOfn[,6])/nsim
    cov.mOfn.50[i] <- sum(psi[i] >= CI.mOfn[,7] & psi[i] <= CI.mOfn[,8])/nsim
    cov.mOfn.25[i] <- sum(psi[i] >= CI.mOfn[,9] & psi[i] <= CI.mOfn[,10])/nsim
    cov.horowitz[i] <- sum(psi[i] >= CI.horowitz[,1] & psi[i] <= CI.horowitz[,2])/nsim
  }
  
  plot(0,0, type = "n", xlim = xregion, ylim = yregion, xaxt = "n", yaxt = "n", xlab = NA, ylab = NA, yaxs = "i", bty = "n")
  rect(lb,0,ub,1, col = "grey90", border = NA)
  par(new = TRUE)
  plot(psi,cov.mOfn.25, xlim = xregion, ylim = yregion, main = NA, type = "l",lty = 4, xlab = NA, ylab = NA, lwd = 2, cex.lab=1.5, cex.axis=1.5, cex.main=1,yaxs = "i", bty = "l")
  #title(ylab = expression(paste("Probability that the CI contains ",psi)), line=2.5, cex.lab=1.5)
  par(new = TRUE)
  ##plot(psi,cov.mOfn.50, xlim = xregion, ylim = yregion, xaxt = "n", yaxt = "n", main = NA, xlab = NA, ylab = NA, type = "l", lty = 5, lwd = 2, cex.lab=1.5, cex.main=1, yaxs = "i", bty = "l")
  ##par(new = TRUE)
  ##plot(psi,cov.mOfn.75, xlim = c(0,1), ylim = c(0,1), xaxt = "n", yaxt = "n", main = NA, xlab = NA, ylab = NA, col = "purple", type = "l", lwd = 1.5, cex.lab=2, cex.main=2, yaxs = "i")
  ##par(new = TRUE)
  ##plot(psi,cov.mOfn.90, xlim = c(0,1), ylim = c(0,1), xaxt = "n", yaxt = "n", main = NA, xlab = NA, ylab = NA, col = "orange", type = "l", lwd = 1.5, cex.lab=2, cex.main=2, yaxs = "i")
  ##par(new = TRUE)
  plot(psi,cov.mOfn.100, xlim = xregion, ylim = yregion, xaxt = "n", yaxt = "n", main = NA, xlab = NA, ylab = NA, type = "l", lty = 3, lwd = 2, cex.lab=1.5, cex.main=1, yaxs = "i", bty = "l")
  par(new = TRUE)
  plot(psi,cov.our, xlim = xregion, ylim = yregion, xaxt = "n", yaxt = "n", main = NA, xlab = NA, ylab = NA, lty = 1, type = "l", lwd = 1.5, cex.lab=1.5, cex.main=1, yaxs = "i", bty = "l")
  par(new = TRUE)
  plot(psi,cov.horowitz, xlim = xregion, ylim = yregion, xaxt = "n", yaxt = "n", main = NA, xlab = NA, ylab = NA, type = "l", lty = 5, lwd = 1.5, cex.lab=1.5, cex.main=1, yaxs = "i", bty = "l")
  abline(h = 0.95, lty = 2)
}


##################################################################################
##MAIN PAPER PLOTS##
##################################################################################

n <- 500

##Setting A (50-50 no restrictions)##
setwd("~/Desktop/plotsForPaper")
pdf("SettingA_N500.pdf")
lb <- 0
ub <- 0.5
nsim <- 5000
homWd <- "~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/binary_noRes_5050"
plotTitle <- "Setting A"
plotCov_singleSetting(n, lb, ub, nsim, homWd, plotTitle, c(0,1), c(0,1))
dev.off()

setwd("~/Desktop/plotsForPaper")
pdf("SettingA_N500_zoomInLeft.pdf")
plotCov_singleSetting_zoomIn(n, lb, ub, nsim, homWd, plotTitle, c(0,0.04), c(0.9,1))
dev.off()

setwd("~/Desktop/plotsForPaper")
pdf("SettingA_N500_zoomInRight.pdf")
plotCov_singleSetting_zoomIn(n, lb, ub, nsim, homWd, plotTitle, c(0.49,0.53), c(0.9,1))
dev.off()

##Setting B (50-50 no harm)##
setwd("~/Desktop/plotsForPaper")
pdf("SettingB_N500.pdf")
lb <- -0.002 ##lb and ub are both 0, but we do -0.002/0.002 so grey rectangle is visible
ub <- 0.002 
nsim <- 5000
homWd <- "~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/binary_noHarm_5050"
plotTitle <- "Setting B"
plotCov_singleSetting(n, lb, ub, nsim, homWd, plotTitle, c(0,1), c(0,1))
dev.off()  

setwd("~/Desktop/plotsForPaper")
pdf("SettingB_N500_zoomIn.pdf")
lb <- -0.0005
ub <- 0.0005
plotCov_singleSetting_zoomIn(n, lb, ub, nsim, homWd, plotTitle, c(0,0.1), c(0.9,1))
dev.off()  


##################################################################################
##SUPP MAT PLOTS##
##################################################################################


##SETTING A##
n <- 200
setwd("~/Desktop/plotsForPaper/suppMatPlots")
pdf("SettingA_N200.pdf")
lb <- 0
ub <- 0.5
nsim <- 5000
homWd <- "~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/binary_noRes_5050"
plotTitle <- "Setting A"
plotCov_singleSetting(n, lb, ub, nsim, homWd, plotTitle, c(0,1), c(0,1))
dev.off()

n <- 1000
setwd("~/Desktop/plotsForPaper/suppMatPlots")
pdf("SettingA_N1000.pdf")
lb <- 0
ub <- 0.5
nsim <- 5000
homWd <- "~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/binary_noRes_5050"
plotTitle <- "Setting A"
plotCov_singleSetting(n, lb, ub, nsim, homWd, plotTitle, c(0,1), c(0,1))
dev.off()

n <- 2000
setwd("~/Desktop/plotsForPaper/suppMatPlots")
pdf("SettingA_N2000.pdf")
lb <- 0
ub <- 0.5
nsim <- 5000
homWd <- "~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/binary_noRes_5050"
plotTitle <- "Setting A"
plotCov_singleSetting(n, lb, ub, nsim, homWd, plotTitle, c(0,1), c(0,1))
dev.off()


##SETTING B##
n <- 200
setwd("~/Desktop/plotsForPaper/suppMatPlots")
pdf("SettingB_N200.pdf")
lb <- -0.002 ##lb and ub are both 0, but we do -0.002/0.002 so grey rectangle is visible
ub <- 0.002 
nsim <- 5000
homWd <- "~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/binary_noHarm_5050"
plotTitle <- "Setting B"
plotCov_singleSetting(n, lb, ub, nsim, homWd, plotTitle, c(0,1), c(0,1))
dev.off()

n <- 1000
setwd("~/Desktop/plotsForPaper/suppMatPlots")
pdf("SettingB_N1000.pdf")
lb <- -0.002 ##lb and ub are both 0, but we do -0.002/0.002 so grey rectangle is visible
ub <- 0.002 
nsim <- 5000
homWd <- "~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/binary_noHarm_5050"
plotTitle <- "Setting B"
plotCov_singleSetting(n, lb, ub, nsim, homWd, plotTitle, c(0,1), c(0,1))
dev.off()

n <- 2000
setwd("~/Desktop/plotsForPaper/suppMatPlots")
pdf("SettingB_N2000.pdf")
lb <- -0.002 ##lb and ub are both 0, but we do -0.002/0.002 so grey rectangle is visible
ub <- 0.002 
nsim <- 5000
homWd <- "~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/binary_noHarm_5050"
plotTitle <- "Setting B"
plotCov_singleSetting(n, lb, ub, nsim, homWd, plotTitle, c(0,1), c(0,1))
dev.off()


##SETTING C##
lb <- 0.25
ub <- 0.5
nsim <- 5000

n <- 200
setwd("~/Desktop/plotsForPaper/suppMatPlots")
pdf("SettingC_N200.pdf")
homWd <- "~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/binary_noRes_7550"
plotTitle <- "Setting C"
plotCov_singleSetting(n, lb, ub, nsim, homWd, plotTitle, c(0,1), c(0,1))
dev.off()

n <- 500
setwd("~/Desktop/plotsForPaper/suppMatPlots")
pdf("SettingC_N500.pdf")
homWd <- "~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/binary_noRes_7550"
plotTitle <- "Setting C"
plotCov_singleSetting(n, lb, ub, nsim, homWd, plotTitle, c(0,1), c(0,1))
dev.off()

n <- 1000
setwd("~/Desktop/plotsForPaper/suppMatPlots")
pdf("SettingC_N1000.pdf")
homWd <- "~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/binary_noRes_7550"
plotTitle <- "Setting C"
plotCov_singleSetting(n, lb, ub, nsim, homWd, plotTitle, c(0,1), c(0,1))
dev.off()

n <- 2000
setwd("~/Desktop/plotsForPaper/suppMatPlots")
pdf("SettingC_N2000.pdf")
homWd <- "~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/binary_noRes_7550"
plotTitle <- "Setting C"
plotCov_singleSetting(n, lb, ub, nsim, homWd, plotTitle, c(0,1), c(0,1))
dev.off()

##SETTING D##
lb <- 0.82
ub <- 0.96
nsim <- 1000

n <- 200
setwd("~/Desktop/plotsForPaper/suppMatPlots")
pdf("SettingD_N200.pdf")
homWd <- "~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/MISTIE_RICV5"
plotTitle <- "Setting D"
plotCov_singleSetting(n, lb, ub, nsim, homWd, plotTitle, c(0,1), c(0,1))
dev.off()

n <- 500
setwd("~/Desktop/plotsForPaper/suppMatPlots")
pdf("SettingD_N500.pdf")
homWd <- "~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/MISTIE_RICV5"
plotTitle <- "Setting D"
plotCov_singleSetting(n, lb, ub, nsim, homWd, plotTitle, c(0,1), c(0,1))
dev.off()

n <- 1000
setwd("~/Desktop/plotsForPaper/suppMatPlots")
pdf("SettingD_N1000.pdf")
homWd <- "~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/MISTIE_RICV5"
plotTitle <- "Setting D"
plotCov_singleSetting(n, lb, ub, nsim, homWd, plotTitle, c(0,1), c(0,1))
dev.off()

n <- 2000
setwd("~/Desktop/plotsForPaper/suppMatPlots")
pdf("SettingD_N2000.pdf")
homWd <- "~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/MISTIE_RICV5"
plotTitle <- "Setting D"
plotCov_singleSetting(n, lb, ub, nsim, homWd, plotTitle, c(0,1), c(0,1))
dev.off()


##################################################################################
##SUPP MAT TABLE OF AVERAGE WIDTH##
##################################################################################

avgWidth <- function(n, nsim, homWd, plotTitle){
  setwd(homWd)
  setwd("ourMethod/results")
  load(paste("res",n,".Rdata", sep = ""))
  if (plotTitle == "Setting D"){
    CI.our <- ci ##CI has previously been computed
  } else {
    CI.our <- matrix(NA, nrow = nsim, ncol = 2) ##compute CI from the confidence set
    for (sim in 1:nsim){
      temp <- which(ci[sim,]==1)
      if (length(temp) >= 1){
        CI.our[sim,] <- c(psi[min(temp)], psi[max(temp)])
      } else{
        CI.our[sim,] <- c(NA,NA)
      }  
    }
  }
  
  setwd(homWd)
  setwd("mOfn/results")
  load(paste("res",n,".Rdata", sep = ""))
  CI.mOfn <- ci[1:nsim,]
  rm(ci)
  
  setwd(homWd)
  setwd("manski_horowitz/results")
  load(paste("res",n,".Rdata", sep = ""))
  CI.horowitz <- ci[1:nsim,]
  rm(ci)
  
  widths <- c(mean(CI.our[,2]-CI.our[,1]),
              mean(CI.mOfn[,2]-CI.mOfn[,1]),
              mean(CI.mOfn[,4]-CI.mOfn[,3]),
              mean(CI.mOfn[,6]-CI.mOfn[,5]),
              mean(CI.mOfn[,8]-CI.mOfn[,7]),
              mean(CI.mOfn[,10]-CI.mOfn[,9]),
              mean(CI.horowitz[,2]-CI.horowitz[,1]))
  return(widths)
}



library(xtable)

##SETTING A##
nsim <- 5000
homWd <- "~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/binary_noRes_5050"
plotTitle <- "Setting A"

n <- 200
settingA.width <- c(n, avgWidth(n, nsim, homWd, plotTitle))
n <- 500
settingA.width <- rbind(settingA.width, c(n, avgWidth(n, nsim, homWd, plotTitle)))
n <- 1000
settingA.width <- rbind(settingA.width, c(n, avgWidth(n, nsim, homWd, plotTitle)))
n <- 2000
settingA.width <- rbind(settingA.width, c(n, avgWidth(n, nsim, homWd, plotTitle)))

##SETTING B##
nsim <- 5000
homWd <- "~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/binary_noHarm_5050"
plotTitle <- "Setting B"

n <- 200
settingB.width <- c(n, avgWidth(n, nsim, homWd, plotTitle))
n <- 500
settingB.width <- rbind(settingB.width, c(n, avgWidth(n, nsim, homWd, plotTitle)))
n <- 1000
settingB.width <- rbind(settingB.width, c(n, avgWidth(n, nsim, homWd, plotTitle)))
n <- 2000
settingB.width <- rbind(settingB.width, c(n, avgWidth(n, nsim, homWd, plotTitle)))


##SETTING C##
nsim <- 5000
homWd <- "~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/binary_noRes_7550"
plotTitle <- "Setting C"

n <- 200
settingC.width <- c(n, avgWidth(n, nsim, homWd, plotTitle))
n <- 500
settingC.width <- rbind(settingC.width, c(n, avgWidth(n, nsim, homWd, plotTitle)))
n <- 1000
settingC.width <- rbind(settingC.width, c(n, avgWidth(n, nsim, homWd, plotTitle)))
n <- 2000
settingC.width <- rbind(settingC.width, c(n, avgWidth(n, nsim, homWd, plotTitle)))


##SETTING D##
nsim <- 1000
homWd <- "~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/MISTIE_RICV5"
plotTitle <- "Setting D"

n <- 200
settingD.width <- c(n, avgWidth(n, nsim, homWd, plotTitle))
n <- 500
settingD.width <- rbind(settingD.width, c(n, avgWidth(n, nsim, homWd, plotTitle)))
n <- 1000
settingD.width <- rbind(settingD.width, c(n, avgWidth(n, nsim, homWd, plotTitle)))
n <- 2000
settingD.width <- rbind(settingD.width, c(n, avgWidth(n, nsim, homWd, plotTitle)))

settingA.width <- data.frame(settingA.width)
settingB.width <- data.frame(settingB.width)
settingC.width <- data.frame(settingC.width)
settingD.width <- data.frame(settingD.width)

lab <- c("n", "our method", "m100", "m90",
         "m75", "m50", "m25", "H-M")
names(settingA.width) <- lab
names(settingB.width) <- lab
names(settingC.width) <- lab
names(settingD.width) <- lab

settingA <- xtable(settingA.width)
digits(settingA) <- c(0,0,rep(2,7))
print(settingA, include.rownames = FALSE)

settingB <- xtable(settingB.width)
digits(settingB) <- c(0,0,rep(2,7))
print(settingB, include.rownames = FALSE)

settingC <- xtable(settingC.width)
digits(settingC) <- c(0,0,rep(2,7))
print(settingC, include.rownames = FALSE)

settingD <- xtable(settingD.width)
digits(settingD) <- c(0,0,rep(2,7))
print(settingD, include.rownames = FALSE)


############################################################################################
##Compute the proportion of trials for which our method yields [0,1] CI##
############################################################################################

span01 <- function(homWd, nsim){
  prop01_ourMethod <- matrix(0, nrow = 4, ncol = 3)
  n <- c(200, 500, 1000, 2000)
  for (i in 1:length(n)){
    load(paste(homWd, "/ourMethod/results/res", n[i], ".Rdata", sep = ""))
    if (setting == "Setting D"){
      CI.our <- ci ##CI has previously been computed
    } else {
      CI.our <- matrix(NA, nrow = nsim, ncol = 2) ##compute CI from the confidence set
      for (sim in 1:nsim){
        temp <- which(ci[sim,]==1)
        if (length(temp) >= 1){
          CI.our[sim,] <- c(psi[min(temp)], psi[max(temp)])
        } else{
          CI.our[sim,] <- c(NA,NA)
        }  
      }
    }
    prop01_ourMethod[i,1] <- max(CI.our[,2]-CI.our[,1]) ##max CI length
    prop01_ourMethod[i,2] <- mean(CI.our[,2]-CI.our[,1]) ##average CI length
    prop01_ourMethod[i,3] <- mean(CI.our[,1] <= 0 & CI.our[,2] >= 1) ##prop of CI's that spanned [0,1]
  }
  
  prop01_mOfn <- matrix(0, nrow = 4, ncol = 5)
  for (i in 1:length(n)){
    load(paste(homWd, "/mOfn/results/res",n[i],".Rdata", sep = ""))
    CI.mOfn <- ci[1:nsim,]  
    
    prop01_mOfn[i,] <- c(max(CI.mOfn[,2]-CI.mOfn[,1]), ##max CI length for the various choices of m
                         max(CI.mOfn[,4]-CI.mOfn[,3]),
                         max(CI.mOfn[,6]-CI.mOfn[,5]),
                         max(CI.mOfn[,8]-CI.mOfn[,7]),
                         max(CI.mOfn[,10]-CI.mOfn[,9]))
    
  }
  
  
  return(list(prop01_ourMethod, prop01_mOfn))
}

##Setting A
setting <- "Setting A"
nsim <- 5000
homWd <- "~/Dropbox/research/github/CI_for_propWhoBenefit/simulations/binary_noRes_5050"
span01(homWd, nsim)

##Setting B
setting <- "Setting B"
nsim <- 5000
homWd <- "~/Dropbox/research/github/CI_for_propWhoBenefit/simulations/binary_noHarm_5050"
span01(homWd, nsim)

##Setting C
setting <- "Setting C"
nsim <- 5000
homWd <- "~/Dropbox/research/github/CI_for_propWhoBenefit/simulations/binary_noRes_7550"
span01(homWd, nsim)

##Setting D
setting <- "Setting D"
nsim <- 1000
homWd <- "~/Dropbox/research/github/CI_for_propWhoBenefit/simulations/MISTIE_RICV5"
span01(homWd, nsim)



###########################################################################################
##Generate the series of boxplots shown in the Supplementary Materials##
###########################################################################################

##These boxplots display the distribution of the observed difference in proportions
##under treatment versus under control, stratified by the left endpoint of the CI
##outputted by our method. This is under Setting B and sample size n = 500.

setwd("~/Dropbox/research/github/CI-for-fraction-who-benefit/simulations/binary_noHarm_5050/ourMethod/results")

load("res500_detailed.Rdata")

result$empdiff <- result$pt2 - result$pc2

boxplot(empdiff ~ ci_left, result, 
        xlab = "Left Endpoint of the Confidence Interval Generated by our Method",
        ylab = "Observed Difference in Proportions")

for (psi in seq(from = 0, to = 0.1, by = 0.01)){
  print(c(psi, 100*mean(result$ci_left == psi)))  
}
