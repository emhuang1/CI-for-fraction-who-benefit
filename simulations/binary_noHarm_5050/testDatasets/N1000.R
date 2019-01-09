rm(list=ls())

nsim <- 5000
set.seed(6881697)
seed <- sample(10^7,nsim)
N <- 1000

datasets <- matrix(NA, nrow = nsim, ncol = 6) ##columns will be pC, pT, pC1, pC2, pT1, pT2
for (sim in 1:nsim){
  set.seed(seed[sim])  
  datasamp <- data.frame(trt = sample(0:1, size = N, replace = TRUE))
  datasamp$y = NA
  NT <- sum(datasamp$trt==1)
  NC <- N - NT
  datasamp$y[datasamp$trt==1] <- sample(1:2, size = NT, replace = TRUE, prob = c(0.5, 0.5))
  datasamp$y[datasamp$trt==0] <- sample(1:2, size = NC, replace = TRUE, prob = c(0.5, 0.5))
  datasets[sim,] <- c(sum(datasamp$trt == 0), sum(datasamp$trt == 1),
                      sum(datasamp$trt == 0 & datasamp$y == 1), sum(datasamp$trt == 0 & datasamp$y == 2),
                      sum(datasamp$trt == 1 & datasamp$y == 1), sum(datasamp$trt == 1 & datasamp$y == 2))
}

datasets <- datasets/N
setwd("~/Dropbox/research/github/CI_for_propWhoBenefit/binary_noHarm_5050/testDatasets")
write.table(datasets, file = "N1000.txt",row.names = FALSE, col.names = FALSE)