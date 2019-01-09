rm(list=ls())

nsim <- 5000
set.seed(6538310)
seed <- sample(10^7,nsim)
N <- 1000

##marg_YC and marg_YT (used below) are vectors of the marginal distributions of RICV5 
##under control and treatment, observed in MISTIE II.

##marg_YC and marg_YT are not shown since the MISTIE II data is proprietary

datasets <- matrix(NA, nrow = nsim, ncol = 14) ##columns will be pC, pT, pC1-pC6, pT1-pT6
for (sim in 1:nsim){
  set.seed(seed[sim])  
  datasamp <- data.frame(trt = sample(0:1, size = N, replace = TRUE))
  datasamp$y = NA
  NT <- sum(datasamp$trt==1)
  NC <- N - NT
  datasamp$y[datasamp$trt==1] <- sample(1:6, size = NT, replace = TRUE, prob = marg_YT)
  datasamp$y[datasamp$trt==0] <- sample(1:6, size = NC, replace = TRUE, prob = marg_YC)
  datasets[sim,] <- c(sum(datasamp$trt == 0), sum(datasamp$trt == 1),
                      sum(datasamp$trt == 0 & datasamp$y == 1), sum(datasamp$trt == 0 & datasamp$y == 2),
                      sum(datasamp$trt == 0 & datasamp$y == 3), sum(datasamp$trt == 0 & datasamp$y == 4),
                      sum(datasamp$trt == 0 & datasamp$y == 5), sum(datasamp$trt == 0 & datasamp$y == 6),
                      sum(datasamp$trt == 1 & datasamp$y == 1), sum(datasamp$trt == 1 & datasamp$y == 2),
                      sum(datasamp$trt == 1 & datasamp$y == 3), sum(datasamp$trt == 1 & datasamp$y == 4),
                      sum(datasamp$trt == 1 & datasamp$y == 5), sum(datasamp$trt == 1 & datasamp$y == 6))
}


setwd("~/Dropbox/research/github/CI_for_propWhoBenefit/MISTIE_RICV5/testDatasets")
write.table(datasets, file = "N1000.txt", row.names = FALSE, col.names = FALSE)
