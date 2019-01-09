rm(list=ls())

setwd("~/Dropbox/research/github/CI_for_propWhoBenefit/CLEARIII/applyHorowitzManski")

#############################################################################
##Libraries##
#############################################################################
library(lpSolveAPI)

#############################################################################
##Function for computing bound estimates##
#############################################################################

boundsNoCov_res <- function(ordinalScale, YT, YC, maxBen, maxHarm){
  nT <- length(YT) #number of treatment subjects
  nC <- length(YC) #number of control subjects
  
  if(nT == 0 | nC == 0){
    #stop("WARNING: YT or YC is empty")
    return(c(NA,NA,NA))
  } else {
    ordinalScale <- sort(ordinalScale, decreasing = FALSE)
    L <- length(ordinalScale)
    varCount <- L^2 #number of pi_{i,j}'s
    
    #The matrix of pi_{i,j}'s is a L x L matrix with varCount pi_{i,j}'s.
    
    scale <- nT * nC ##this can be made large to help with solving the linear program
    
    #############################################################################
    ##Calculate marginal cdf's##
    #############################################################################
    cdf_C <- rep(0,L)
    for (i in 1:L){
      cdf_C[i] <- sum(YC <= ordinalScale[i])/nC
    }
    cdf_T <- rep(0,L)
    for (i in 1:L){
      cdf_T[i] <- sum(YT <= ordinalScale[i])/nT
    }
    
    #############################################################################
    ##Which pi_{i,j}'s are affected by the restrictions?##
    #############################################################################
    restrictions <- matrix(0,nrow=L,ncol=L)
    for (i in 1:L){
      for (j in 1:L){
        if (ordinalScale[j]-ordinalScale[i]>maxBen | ordinalScale[i]-ordinalScale[j]>maxHarm){
          restrictions[i,j] <- 1
        }
      }
    }
    restrictions <- c(restrictions)
    #############################################################################
    ##LP for epsilon##
    #############################################################################
    lprec <- make.lp(0,(varCount+1))
    
    ##Setting objective function
    objfn <- c(rep(0,varCount),1)
    set.objfn(lprec,objfn)
    
    ##Setting non-negativity bounds
    set.bounds(lprec, lower = rep(0, (varCount+1)), upper = NULL)
    
    ##pi_{i,j}'s sum to 1
    add.constraint(lprec, xt = c(rep(1,varCount),0), "=", rhs = scale)
    
    ##incorporating the restrictions
    add.constraint(lprec, xt = c(restrictions,0), "=", rhs = 0)
    
    ##marginal cdf constraints
    pij.matrix <- matrix(1:varCount, nrow = L, ncol = L)
    
    for (i in 1:(L-1)){
      add.constraint(lprec, xt = c(rep(1,L*i),-1), "<=", rhs = (cdf_C[i]*scale), indices = c(pij.matrix[1:i,],(varCount+1))) 
      add.constraint(lprec, xt = rep(1,L*i+1), ">=", rhs = (cdf_C[i]*scale), indices = c(pij.matrix[1:i,],(varCount+1)))
    }
    
    for (i in 1:(L-1)){
      add.constraint(lprec, xt = c(rep(1,L*i),-1), "<=", rhs = (cdf_T[i]*scale), indices = c(pij.matrix[,1:i],(varCount+1))) 
      add.constraint(lprec, xt = rep(1,L*i+1), ">=", rhs = (cdf_T[i]*scale), indices = c(pij.matrix[,1:i],(varCount+1)))
    }
    
    ##Solving linear program
    eps.flag <- solve(lprec)
    
    if(eps.flag != 0){
      stop("WARNING: problem with LP for eps")
    }
    
    eps <- get.objective(lprec)/scale
    
    if (eps < 10 ^ (-10)){
      eps <- 0
    }
    
    rm(lprec)
    
    #############################################################################
    ##LP for lb and ub##
    #############################################################################
    
    lprec <- make.lp(0,varCount)
    
    ##Setting objective function
    objfn <- matrix(0, L, L)
    for (i in 1:L){
      for (j in 1:L){
        if (ordinalScale[j]>ordinalScale[i]){ 
          objfn[i,j] <- 1
        }
      }
    }
    objfn <- c(objfn)
    set.objfn(lprec,objfn)
    
    ##Setting non-negativity bounds 
    set.bounds(lprec, lower = rep(0, varCount), upper = NULL) 
    
    ##pi_{i,j}'s sum to 1
    add.constraint(lprec, xt = rep(1,varCount), "=", rhs = scale)
    
    ##incorporating the restrictions
    add.constraint(lprec, xt = restrictions, "=", rhs = 0)
    
    ##marginal cdf constraints
    for (i in 1:(L-1)){
      add.constraint(lprec, xt = rep(1,L*i), "<=", rhs = (cdf_C[i]+eps)*scale, indices = c(pij.matrix[1:i,])) 
      add.constraint(lprec, xt = rep(1,L*i), ">=", rhs = (cdf_C[i]-eps)*scale, indices = c(pij.matrix[1:i,]))
    }
    
    for (i in 1:(L-1)){
      add.constraint(lprec, xt = rep(1,L*i), "<=", rhs = (cdf_T[i]+eps)*scale, indices = c(pij.matrix[,1:i])) 
      add.constraint(lprec, xt = rep(1,L*i), ">=", rhs = (cdf_T[i]-eps)*scale, indices = c(pij.matrix[,1:i]))
    }
    
    lb.flag <- solve(lprec)
    
    if(lb.flag != 0){
      stop("WARNING: problem with LP for lb")
    }
    
    lb <- get.objective(lprec)/scale
    
    lp.control(lprec,sense='max')
    
    ub.flag <- solve(lprec)
    
    if(ub.flag != 0){
      stop("WARNING: problem with LP for ub")
    }
    
    ub <- get.objective(lprec)/scale
    
    return(c(lb,ub,eps))
  }
}

###########################################################################
##Function to Construct 95% CI using method by Horowitz and Manski (2000)##
###########################################################################

##maxBen, maxHarm are restrictions on how much a person can benefit or be harmed
##by treatment, compared to control

##outcomeLevels are the possible values of the outcome

##datasamp is the data set with "trt" column (0 = control, 1 = treatment) and outcome "y" column 

getCI.HM <- function(maxBen, maxHarm, outcomeLevels, datasamp, B){
  ##Bound estimates for trial data
  bounds_trial <- boundsNoCov_res(outcomeLevels, datasamp$y[datasamp$trt==1], datasamp$y[datasamp$trt==0], maxBen, maxHarm)
  lb_trial <- bounds_trial[1]
  ub_trial <- bounds_trial[2]
  
  ##Sample size
  N <- nrow(datasamp)
  
  ##Bound estimates for bootstrap replicate data sets
  bounds <- matrix(NA, nrow = B, ncol = 3)
  for (b in 1:B){
    dataresamp <- datasamp[sample(1:N,N, replace = TRUE),]
    YT <- dataresamp$y[dataresamp$trt == 1]
    YC <- dataresamp$y[dataresamp$trt == 0]
    bounds[b,] <- boundsNoCov_res(outcomeLevels, YT, YC, maxBen, maxHarm)
  }
  
  lb_star <- bounds[,1]
  ub_star <- bounds[,2]
  
  ##Find z_na
  z_candidates <- seq(from = 0, to = 1, by = 0.001) ##candidates for z_na
  prob_z <- sapply(X = z_candidates, FUN = function(z){mean(lb_star - z <= lb_trial & ub_trial <= ub_star + z)})
  z_na <- min(z_candidates[which(prob_z >= 0.95)])
  
  ##Get CI
  leftLimit <- lb_trial - z_na
  rightLimit <- ub_trial + z_na
  
  if (leftLimit < 0){ ##We know the lower bound is 0 at minimum
    leftLimit <- 0
  }
  
  if (rightLimit > 1){ ##We know the upper bound is 1 at maximum
    rightLimit <- 1
  }
  
  CI <- c(leftLimit, rightLimit)
  return(CI)
}


##########################################################################
##Load Data and Set Inputs##
##########################################################################

setwd("~/Dropbox/research/github/CI_for_propWhoBenefit/CLEARIII/constructDataset")
load("data.Rdata")

##No restrictions
maxBen <- 100
maxHarm <- 100

##Number of bootstrap replicate data sets
B <- 10000

##########################################################################
##30-day mortality##
##########################################################################
mort30data <- data.frame(trt = data$tmt, y = -data$mortality30)
mort30.ci <- getCI.HM(maxBen, maxHarm, c(-1,0), mort30data, B)

##########################################################################
##180-day mortality##
##########################################################################
mort180data <- data.frame(trt = data$tmt, y = -data$mortality180)
mort180data <- subset(mort180data, !is.na(y))
mort180.ci <- getCI.HM(maxBen, maxHarm, c(-1,0), mort180data, B)

##########################################################################
##30-day mRS##
##########################################################################
mrs30data <- data.frame(trt = data$tmt, y = -data$rankin30)
mrs30data <- subset(mrs30data, !is.na(y))
mrs30.ci <- getCI.HM(maxBen, maxHarm, (-6):0, mrs30data, B)

##########################################################################
##180-day mRS##
##########################################################################
mrs180data <- data.frame(trt = data$tmt, y = -data$rankin180)
mrs180data <- subset(mrs180data, !is.na(y))
mrs180.ci <- getCI.HM(maxBen, maxHarm, (-6):0, mrs180data, B)

round(mrs30.ci, digits = 2)
round(mrs180.ci, digits = 2)
round(mort30.ci, digits = 2)
round(mort180.ci, digits = 2)
