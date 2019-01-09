rm(list=ls())

setwd("~/Dropbox/research/github/CI_for_propWhoBenefit/CLEARIII/applyMoutOfN")

source("mOfn_code.R")

maxBen <- 100
maxHarm <- 100


setwd("~/Dropbox/research/github/CI_for_propWhoBenefit/CLEARIII/constructDataset")
load("data.Rdata")

##########################################################################
##30-day mortality##
##########################################################################
mort30data <- data.frame(trt = data$tmt, y = -data$mortality30)
mort30.ci <- getCI.mOfn(maxBen, maxHarm, c(-1,0), mort30data)

##########################################################################
##180-day mortality##
##########################################################################
mort180data <- data.frame(trt = data$tmt, y = -data$mortality180)
mort180data <- subset(mort180data, !is.na(y))
mort180.ci <- getCI.mOfn(maxBen, maxHarm, c(-1,0), mort180data)

##########################################################################
##30-day mRS##
##########################################################################
mrs30data <- data.frame(trt = data$tmt, y = -data$rankin30)
mrs30data <- subset(mrs30data, !is.na(y))
mrs30.ci <- getCI.mOfn(maxBen, maxHarm, (-6):0, mrs30data)

##########################################################################
##180-day mRS##
##########################################################################
mrs180data <- data.frame(trt = data$tmt, y = -data$rankin180)
mrs180data <- subset(mrs180data, !is.na(y))
mrs180.ci <- getCI.mOfn(maxBen, maxHarm, (-6):0, mrs180data)

round(mrs30.ci, digits = 2)
round(mrs180.ci, digits = 2)
round(mort30.ci, digits = 2)
round(mort180.ci, digits = 2)
