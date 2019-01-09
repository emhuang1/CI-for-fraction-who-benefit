##6/24/17

##Plot distribution of T

library(R.matlab)

rm(list=ls())


plotResult <- function(data){
  tdraws <- data$tdraws
  psi <- as.numeric(data$psi)
  
  for (i in 1:length(psi)){
    hist(tdraws[i,], breaks = 50, main = psi[i], probability = TRUE)
  }
}


###########################################################################
##Figure 2 in paper; Setting A, psi = 0.5##

setwd("~/Dropbox/research/github/CI_for_propWhoBenefit/limitDist_underNull/binary_noRes_5050")

data <- readMat("result23-Jun-2017.mat")

tdraws <- data$tdraws
psi <- as.numeric(data$psi)

setwd("~/Dropbox/research/github/CI_for_propWhoBenefit/limitDist_underNull")
pdf("Fig2.pdf")
hist(tdraws[5,], probability = TRUE, breaks = 100, xlab = "", main = "")
dev.off()

mean(tdraws[5,]<= 10^(-10))
##0.25105
