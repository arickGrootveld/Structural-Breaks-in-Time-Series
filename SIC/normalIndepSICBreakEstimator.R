############################ normalIndepSICBreakEstimator.R ##################################
# Function utilizing SIC calculations to detect statistical breaks
#
# v1.0.0
# Contributors: Jamie Shannon, Arick Grootveld
#############################################################################

## Supporting functions used by the main function

# From Noah Prime's paper the 1st summation term (page 6 middle of the page)
minsum<-function(z,y,d){
  sum(((z-d)^2)[c(1:y)])
}
# From Noah Prime's paper the 2nd summation term (page 6 middle of the page)
minsum2<-function(h,j,k,l){
  sum(((h-j)^2)[c(k:l)])
}
# Calculating the null hypothesis value
SICNull<-function(a,b) {
  a*log(2*pi)+b+log(a)
}
# Calculating the alternate hypothesis value
SICalt<-function(c,d) {
  c*log(2*pi)+d+2*log(c)
}

normalIndepSICCalc <- function(simMatrix, param1, param2){
  simDims <- dim(simMatrix)
  Cor <- 0
  detectedBreakIndexes = matrix(-1, 1, simDims[1])
  # Looping through each row of the simulation matrix, getting a single simulation to evaluate our performance on
  for (simNum in 1:simDims[1]){
    Dist <- simMatrix[simNum,] # The target data to detect a break from
    n = length(Dist) # Number of samples in this instance of the data
    Signull=sum((Dist-mean(Dist))^2)
    g1 <- -1 # Initializing the SIC test statistic
    # Looping through each index so we can get a simple 
    for (k in 1:(n-1)){
      g<-(minsum(Dist,k,param1)+minsum2(Dist,param2,(k+1),n))
      if(k == 1){
        g1 <- g
        potentialIndex = k-1
      } else if (g<= g1){
        g1=g
        potentialIndex = k-1
      }
    }
    # Comparison between the null and alternative
    if((SICalt(n,g1)) < SICNull(n,Signull)){
      Cor=Cor+1
      detectedBreakIndexes[1, simNum] = potentialIndex # Since a change was detected we grab the index we detected it at
    }
  }
  coverageProb <- Cor/simDims[1]
  return(c(coverageProb, detectedBreakIndexes))
}