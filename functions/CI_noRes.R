rm(list=ls())

################################################################################################
##Libraries##
################################################################################################

library(mvtnorm)
library(Matrix)
library(gurobi)

################################################################################################
##GIVE THESE INPUTS##
################################################################################################

#Below are inputs from an example data set. You can use these inputs to practice using our code. Update the inputs to get the confidence interval for your own data set.

#Sample size
n <- 500

#Number of levels of ordinal outcome
L <- 2

#Row vector of number of participants assigned to control and with outcome = 1, number of participants assigned to control and with outcome = 2, ...., number of participants assigned to control and with outcome = L

nCvec <- c(144,113)

#Row vector of number of participants assigned to treatment and with outcome = 1, number of participants assigned to treatment and with outcome = 2, ...., number of participants assigned to treatment and with outcome = L

nTvec <- c(56, 187)

#Define nu (see Appendix J of Supplementary Materials)
nu <- (100 * sqrt(n)/log(n))^2

################################################################################################

pC <- sum(nCvec)/n
pT <- sum(nTvec)/n
pCvec <- nCvec/n
pTvec <- nTvec/n

#Options for QP solver
params <- list(OutputFlag = 0, OptimalityTol=1e-9,  BarConvTol=1e-10) 
##does not print output
##sets tolerances

psi <- seq(from = 0, to = 1, by = 0.01) 

numGam <- 2*L #these numbers will be used often
numPI <- L^2
numBoth <- numGam + numPI

j <- (1:numPI)+numGam
Amat1 <- sparseMatrix(i = 1:numPI,j = j, x = -1, dims = c(numPI,numBoth))
bvec1 <- rep(0,numPI)
dvec1 <- -2 * c(pCvec, pTvec, rep(0,numPI))
v = c(rep(pC, L),rep(pT, L))*2
Dmat1 <- sparseMatrix(i = 1:numGam, j = 1:numGam, x = v, dims = c(numBoth,numBoth))

mat <- t(matrix(1:numPI, nrow = L, ncol = L))
benPI <- mat[upper.tri(mat)]
i1 <- rep(1,numPI)
j1 <- (1:numPI) + numGam
i2 <- 2*rep(1,L*(L-1)/2)
j2 <- numGam + benPI
i3 <- rep(1:numGam,each = L) + 2
j3 <- c(1:numPI,as.numeric(mat))+numGam
i4 <- (1:numGam)+2
j4 <- 1:numGam
i <- c(i4,i1,i2,i3)
j <- c(j4,j1,j2,j3)
v <- rep(1,length(i))
v[1:length(i4)] <- -1
Aeq1 <- sparseMatrix(i = i, j = j, x = v, dims = c(2+numGam,numBoth))

rm(i,j,i1,i2,i3,i4,j1,j2,j3,j4) 

bvec2 <- rep(0,numGam)
Amat2 <- sparseMatrix(i = 1:numGam,j = 1:numGam,x = -1,dims = c(numGam,numGam))
beq2 <- t(c(1,1))
Aeq2 <- rbind(c(rep(1,L),rep(0,L)), c(rep(0,L),rep(1,L)))
dvec2 <- -2*c(pCvec, pTvec)
v <- c(rep(pC, L),rep(pT, L))*2
Dmat2 <- sparseMatrix(i = 1:numGam,j = 1:numGam,x = v,dims = c(numGam,numGam))

model       <- c()
model$Q     <- Dmat2/2
model$obj   <- dvec2
model$A     <- rbind(Amat2,Aeq2)
model$sense <- c(rep('<',dim(Amat2)[1]),rep('=',dim(Aeq2)[1]))
model$rhs   <- c(bvec2,beq2)
model$lb    <- -Inf

result   <- gurobi(model, params)
gammaHat <- result$x
val2     <- result$objval

val2 <- val2 + 1


sigma <- matrix(0, ncol = numGam, nrow = numGam)
v <- c(rep(pC, L),rep(pT, L))
j <- c(pCvec, pTvec)
for (r in 1:numGam){
  for (c in 1:numGam){
    if (r == c){
      sigma[r,r] <- j[r] - 2*gammaHat[r]*j[r]+gammaHat[r]^2*v[r]
    } else if (r > L && c <= L){
      sigma[r,c] <- 0  
    } else if (c > L && r <= L){
      sigma[r,c] <- 0
    } else {
      sigma[r,c] <- -gammaHat[r]*j[c]-gammaHat[c]*j[r] + gammaHat[c]*gammaHat[r]*v[r]
    }
  }
}
sigma <- sigma * 4



bvec34 <- c(rep(0,numPI),sqrt(nu))
beq4 <- rep(0,numGam * 2)
beq3 <- rep(0,numGam * 2 + 1)

Amat34 <- sparseMatrix(i = c(1:numPI,rep(numPI+1,numPI)),
                       j = c(1:numPI,1:numPI),
                       x = c(rep(-1,numPI),rep(1,numPI)),
                       dims = c(numPI+1,numPI+numGam*2))
Dmat34 <- sparseMatrix(i = (1:numGam)+numPI,j = (1:numGam)+numPI,x = 1,dims = c(numPI+numGam*2,numPI+numGam*2))

i1 <- rep(1:numGam, each = L)
j1 <- c(1:numPI,as.numeric(mat))
v1 <- rep(1,length(i1))
i2 <- 1:(numGam*2)
j2 <- c((1:numGam)+numPI + numGam,(1:numGam)+numPI + numGam)
v2 <- -1 * rep(1,length(i2))
i3 <- (numGam+1):(2*numGam)
v3 <- rep(1,length(i3))
j3 <- (1:numGam)+numPI
i4 <- rep(i3,numPI)
j4 <- rep(1:numPI,each = numGam)
v4 <- rep(gammaHat,numPI)
i <- c(i1,i2,i3,i4)
j <- c(j1,j2,j3,j4)
v <- c(v1,v2,v3,v4)
Aeq4 <- sparseMatrix(i = i,j = j,x = v,dims = c(numGam*2,numPI+numGam*2))
rm(i1,i2,i3,i4,j1,j2,j3,j4,v1,v2,v3,v4) 

ndraw <- 1000
T95 <- rep(0,length(psi))

##iPsi,jPsi are used to make Aeq3 inside the loop
jPsi <- c(j,1:numPI,benPI)
iPsi <- c(i,(1+numGam*2)*rep(1,numPI+length(benPI)))

rm(i,j)

for (p in 1:length(psi)){
  tdraws <- rep(0,ndraw)
  
  vPsi <- c(v,psi[p]*rep(1,length(1:numPI)),-1*rep(1,length(benPI)))
  Aeq3 <- sparseMatrix(i = iPsi,j = jPsi,x = vPsi,dims = c(1+numGam*2,numPI+numGam*2))
  
  for (d in 1:ndraw){
    Z <- rmvnorm(n = 1, mean = rep(0,numGam), sigma = sigma)
    dvec34 <- c(rep(0,numPI),Z,rep(0,numGam))
    model       <- c()
    model$Q     <- Dmat34/2
    model$obj   <- dvec34
    model$A     <- rbind(Amat34,Aeq4)
    model$sense <- c(rep('<',dim(Amat34)[1]),rep('=',dim(Aeq4)[1]))
    model$rhs   <- c(bvec34,beq4)
    model$lb    <- -Inf
    
    result   <- gurobi(model, params)
    val4     <- result$objval
    
    model       <- c()
    model$Q     <- Dmat34/2
    model$obj   <- dvec34
    model$A     <- rbind(Amat34,Aeq3)
    model$sense <- c(rep('<',dim(Amat34)[1]),rep('=',dim(Aeq3)[1]))
    model$rhs   <- c(bvec34,beq3)
    model$lb    <- -Inf
    
    result   <- gurobi(model, params)
    val3     <- result$objval
    
    tdraws[d] <- val3 - val4
  }

  tdraws <- sort(tdraws)
  T95[p] <- tdraws[0.95*ndraw]
  
  beq1 <- c(1,psi[p], rep(0,numGam))
  model       <- c()
  model$Q     <- Dmat1/2
  model$obj   <- dvec1
  model$A     <- rbind(Amat1,Aeq1)
  model$sense <- c(rep('<',dim(Amat1)[1]),rep('=',dim(Aeq1)[1]))
  model$rhs   <- c(bvec1,beq1)
  model$lb    <- -Inf
  
  result   <- gurobi(model, params)
  gammaHat <- result$x
  val1     <- result$objval
  
  val1 <- val1 + 1
  
  Tn <- n * (val1 - val2)
  
  if (Tn <= T95[p]+1e-10){
    break
  }

}


leftLim <- psi[p]

for (p in seq(from = length(psi), to = 1, by = -1)){
  tdraws <- rep(0,ndraw)
  
  vPsi <- c(v,psi[p]*rep(1,length(1:numPI)),-1*rep(1,length(benPI)))
  Aeq3 <- sparseMatrix(i = iPsi,j = jPsi,x = vPsi,dims = c(1+numGam*2,numPI+numGam*2))
  
  for (d in 1:ndraw){
    Z <- rmvnorm(n = 1, mean = rep(0,numGam), sigma = sigma)
    dvec34 <- c(rep(0,numPI),Z,rep(0,numGam))
    model       <- c()
    model$Q     <- Dmat34/2
    model$obj   <- dvec34
    model$A     <- rbind(Amat34,Aeq4)
    model$sense <- c(rep('<',dim(Amat34)[1]),rep('=',dim(Aeq4)[1]))
    model$rhs   <- c(bvec34,beq4)
    model$lb    <- -Inf
    
    result   <- gurobi(model, params)
    gammaHat <- result$x
    val4     <- result$objval
    
    model       <- c()
    model$Q     <- Dmat34/2
    model$obj   <- dvec34
    model$A     <- rbind(Amat34,Aeq3)
    model$sense <- c(rep('<',dim(Amat34)[1]),rep('=',dim(Aeq3)[1]))
    model$rhs   <- c(bvec34,beq3)
    model$lb    <- -Inf
    
    result   <- gurobi(model, params)
    gammaHat <- result$x
    val3     <- result$objval
    
    tdraws[d] <- val3 - val4
  }
  
  tdraws <- sort(tdraws)
  T95[p] <- tdraws[0.95*ndraw]
  
  beq1 <- c(1,psi[p],rep(0,numGam))
  
  model       <- c()
  model$Q     <- Dmat1/2
  model$obj   <- dvec1
  model$A     <- rbind(Amat1,Aeq1)
  model$sense <- c(rep('<',dim(Amat1)[1]),rep('=',dim(Aeq1)[1]))
  model$rhs   <- c(bvec1,beq1)
  model$lb    <- -Inf
  
  result   <- gurobi(model, params)
  gammaHat <- result$x
  val1     <- result$objval
  val1 <- val1 + 1
  
  Tn <- n * (val1 - val2)
  
  if (Tn <= T95[p]+1e-10){
    break
  }

}

rightLim <- psi[p]

CI <- c(leftLim, rightLim)


################################################################################################
##The confidence interval for the example data set is [0.27,0.60].##
################################################################################################

