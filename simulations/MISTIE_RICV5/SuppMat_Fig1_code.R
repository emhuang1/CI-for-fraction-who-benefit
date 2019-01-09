#######################################################################################
##Function for getting RICV5 from the continuous reduction in clot volume##
#######################################################################################

levelfun <- function(x){ 
  if(x < 0){
    y = 1
  } else if(x >= 0 & x < 5){
    y = 2
  } else if(x >= 5 & x < 10){
    y = 3
  } else if(x >= 10 & x < 15){
    y = 4
  } else if(x >= 15 & x < 20){
    y = 5
  } else {
    y = 6
  } 
  return(y)
}

#######################################################################################
##Compute marginals##
#######################################################################################

load("MISTIEII_deidentified.Rdata") 
##This data set is not included with the paper since it is proprietary
##Below we get the outcome RICV5 for the treatment arm participants (YTd)
##and control arm participants (YCd), separately.

data <- subset(data,Follow_up_Visit == 30)
data <- data.frame(id = data$patientName, tmt = data$Group_Assigned, init = data$Pre_Rand_ICH_Volume_RC, post = data$eot_ich_9_13)
data$volChange <- data$init - data$post
control <- subset(data, tmt == "Medical")
tmt <- subset(data, tmt == "Surgical")
YT <- tmt$volChange
YC <- control$volChange
YCd <- sapply(YC, levelfun)
YTd <- sapply(YT, levelfun)

table(YCd)
table(YTd)
nrow(control)
nrow(tmt)


########################################################################################
##Plot of Marginal Distributions (Figure 1 in Supplementary Materials)##
########################################################################################
nLevels <- 6
countC <- rep(0,nLevels)
countT <- rep(0,nLevels)

for(i in 1:nLevels){
  countC[i] <- sum(YCd==i)
  countT[i] <- sum(YTd==i)
}

nT <- length(YTd)
nC <- length(YCd)

setwd("~/Dropbox/research/github/CI_for_propWhoBenefit/MISTIE_RICV5")

pdf("RICV5.pdf")
mat <- rbind(countT, countC)
colnames(mat) <- 1:nLevels
rownames(mat) <- c("treatment", "control")
prop <- prop.table(mat, margin = 1)
barplot(prop, main="", col=c("red","blue"), beside=TRUE, ylim = c(0,1),yaxt="n",
        xlab = "RICV5", ylab = NA, cex.main = 2, cex.lab = 2, cex.axis = 2, cex.names = 2)
axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1), labels = c("0","0.2","0.4","0.6","0.8","1"),las = 2, cex.axis = 2)
text(x = 5, y = 0.9, labels = bquote(n[T]~"="~.(nT)~","~n[C]~"="~.(nC)),cex=2)
legend(x=12,y=1,legend = rownames(mat), fill = c("red", "blue"),cex = 1.5)
dev.off()