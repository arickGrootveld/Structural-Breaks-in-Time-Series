############################ CUSUMCalc.R ##################################
# Function utilizing CUSUM calculations to detect statistical breaks
#
# v1.0.0
# Contributors: Jamie Shannon, Arick Grootveld
#############################################################################
CUSUMCalc <- function(simMatrix, critVal=0.0133, longRunVar=1){
  simDims = dim(simMatrix)

  Cor=0 #Count for the number of correct changes detected
  Mty=NULL
  # Preallocating an array to contain the indexes of the estimated breaks, so we can calculate the break point location accuracy of
  # our algorithms
  detectedBreakIndexes = matrix(-1, 1, dim(simMatrix)[1])
  
  for(simNum in c(1:simDims[1])){
    Dist <- simMatrix[simNum,] # The target data to detect a break from
    Mn=0 #setting the maxtype stat to 0
    k=1 #Set k to 1 since we are checking all k's in the range {1,...,n}
    Ts=0 #Set the test statistic equal to 0 which we will compare to Mn
    Mt=0
    while(k<=length(Dist)){ #While loop runs while k is in the range {1,...,n}
      # Using absolute value since we only care about the difference between the long run and short run means
      # not the signed difference
      Ts<-abs(TestStat(Dist, k/(length(Dist)), length(Dist))) #Setting the teststat equal to Ts for different k's
      k=k+1 #Adding to k after each iteration to keep the loop going
      if(Ts>=Mn){ #Comparing Ts to Mn
        Mn=Ts #If the Ts is bigger than Mn will be the "max" of all the k's tested so far
        # Setting a variable to grab the index with the best chance of being a break location to store if we detect a break
        potentialIndex = k-1
      }
      Mt=MaxType(Mn,longRunVar) #Setting the Max Type equal to Mt
    }
    if(Mt>critVal){ #Comparing the Max Type to Asymptotic value
      Cor=Cor+1 #Counting this as a change we detected
      detectedBreakIndexes[1, simNum] = potentialIndex # Since a change was detected we grab the index we detected it at
    }
    Mty=c(Mty, MaxType(Mn,longRunVar))
  }
  coverageProbability <- Cor/simDims[1]
  # Just returning the coverage probabilitiy for right now, may want to return other stuff later as well
  return(c(coverageProbability, detectedBreakIndexes))
}

# Moved function declarations outside of loops
# Should probably talk to the group about this
TestStat<-function(y,x,n){ #Test Statistic (Zn) from Aue and ALexander
  return((1/(sqrt(n)))*((sum(y[1:floor(n*x)]))-(((floor(n*x))/(n))*(sum(y[1:n])))))
}

MaxType=function(m,w){ #Max Type Function from Aue Alexander
  return((1/w)*m)
}