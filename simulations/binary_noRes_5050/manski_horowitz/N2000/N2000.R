##10/19/18

##Implement Manski-Horowitz

rm(list=ls())

val <- as.numeric(Sys.getenv("SGE_TASK_ID"))

nsim <- 5000
set.seed(9814663)
seed <- sample(10^7,nsim)


#############################################################################
##Libraries##
#############################################################################
library(lpSolveAPI)
library(boot)

#############################################################################
##Function for estimating bounds##
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



#############################################################################
##Inputs##
#############################################################################
maxBen <- 100
maxHarm <- 100

N <- 2000 ##update

#############################################################################
##Simulate randomized trial##
#############################################################################
B <- 10000
CI <- rep(NA,2)

set.seed(seed[val])
datasamp <- data.frame(trt = sample(0:1, size = N, replace = TRUE))
datasamp$y = NA
NT <- sum(datasamp$trt==1)
NC <- N - NT
datasamp$y[datasamp$trt==1] <- sample(1:2, size = NT, replace = TRUE, prob = c(0.5, 0.5))
datasamp$y[datasamp$trt==0] <- sample(1:2, size = NC, replace = TRUE, prob = c(0.5, 0.5))

#############################################################################
##Get CI's##
#############################################################################

##Bound estimates for trial data
bounds_trial <- boundsNoCov_res(1:2, datasamp$y[datasamp$trt==1], datasamp$y[datasamp$trt==0], maxBen, maxHarm)
lb_trial <- bounds_trial[1]
ub_trial <- bounds_trial[2]

##Bound estimates for bootstrap replicate data sets
bounds <- matrix(NA, nrow = B, ncol = 3)
for (b in 1:B){
  dataresamp <- datasamp[sample(1:N,N, replace = TRUE),]
  YT <- dataresamp$y[dataresamp$trt == 1]
  YC <- dataresamp$y[dataresamp$trt == 0]
  bounds[b,] <- boundsNoCov_res(1:2, YT, YC, maxBen, maxHarm)
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

save(CI, file=paste("N2000","-","seed", val, "-", format(Sys.time(), "%Y%m%d-%H%M"),".Rdata", sep=""))

