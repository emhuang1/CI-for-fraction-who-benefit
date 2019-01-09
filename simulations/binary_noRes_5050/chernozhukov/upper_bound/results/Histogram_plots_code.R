##9/18/18

##For each sample size n, plot a histogram of the right endpoints of the 
##97.5% one-sided CI's for the upper bound

##Check the coverage probability of the upper bound, which is 0.5 in this simulation setting 
##(The coverage probability should be 0.975 or higher)

rm(list=ls())

##Set working directory to where the results files are saved
setwd("~/Dropbox/research/github/CI_for_propWhoBenefit/simulations/binary_noRes_5050/chernozhukov/upper_bound/results")

##Store plots in a pdf file
pdf("Histogram_plots.pdf")

#######################################################################################
##N = 200##
#######################################################################################

load("res200.Rdata")

##Plot histogram
N <- 200
h <- hist(theta, breaks = 50, plot = FALSE) 
h$density <- h$counts/sum(h$counts)
plot(h,freq=FALSE,
     xlab = "Upper Endpoint of One-sided CI for Upper Bound",
     ylab = "Proportion of Simulations",
     xlim = c(0, 1),
     main = paste("Histogram of CI Endpoint for Upper Bound\n", "n = ", N))
abline(v = 0.5, col = "red")

##Compute coverage probability of upper bound
print(mean(theta >= 0.5))

#######################################################################################
##N = 500##
#######################################################################################

load("res500.Rdata")

##Plot histogram
N <- 500
h <- hist(theta, breaks = 50, plot = FALSE) 
h$density <- h$counts/sum(h$counts)
plot(h,freq=FALSE,
     xlab = "Upper Endpoint of One-sided CI for Upper Bound",
     ylab = "Proportion of Simulations",
     xlim = c(0, 1),
     main = paste("Histogram of CI Endpoint for Upper Bound\n", "n = ", N))
abline(v = 0.5, col = "red")

##Compute coverage probability of upper bound
print(mean(theta >= 0.5))

#######################################################################################
##N = 1000##
#######################################################################################

load("res1000.Rdata")

##Plot histogram
N <- 1000
h <- hist(theta, breaks = 50, plot = FALSE) 
h$density <- h$counts/sum(h$counts)
plot(h,freq=FALSE,
     xlab = "Upper Endpoint of One-sided CI for Upper Bound",
     ylab = "Proportion of Simulations",
     xlim = c(0, 1),
     main = paste("Histogram of CI Endpoint for Upper Bound\n", "n = ", N))
abline(v = 0.5, col = "red")

##Plot coverage probability of upper bound
print(mean(theta >= 0.5))

#######################################################################################
##N = 2000##
#######################################################################################

load("res2000.Rdata")

##Plot histogram
N <- 2000
h <- hist(theta, breaks = 50, plot = FALSE) 
h$density <- h$counts/sum(h$counts)
plot(h,freq=FALSE,
     xlab = "Upper Endpoint of One-sided CI for Upper Bound",
     ylab = "Proportion of Simulations",
     xlim = c(0, 1),
     main = paste("Histogram of CI Endpoint for Upper Bound\n", "n = ", N))
abline(v = 0.5, col = "red")

##Compute coverage probability of the upper bound
print(mean(theta >= 0.5))

dev.off()

