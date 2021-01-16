CUSUMCalc <- function(simMatrix, critVal=2.408, longRunVar=1){
  simDims = dim(simMatrix)

  Cor=0 #Count for the number of correct changes detected
  Mty=NULL
  
  for(simNum in c(1:simDims[1])){
    Dist <- simMatrix[simNum,] #A normal distribution where the mean changes from 5 to 7 at the 20th observation
    Mn=0 #setting the maxtype stat to 0
    k=1 #Set k to 1 since we are checking all k's in the range {1,...,n}
    Ts=0 #Set the test statistic equal to 0 which we will compare to Mn
    Mt=0
    while(k<=length(Dist)){ #While loop runs while k is in the range {1,...,n}
      Ts<-TestStat(Dist, k/(length(Dist)), length(Dist)) #Setting the teststat equal to Ts for different k's
      k=k+1 #Adding to k after each iteration to keep the loop going
      if(Ts>=Mn){ #Comparing Ts to Mn
        Mn=Ts #If the Ts is bigger than Mn will be the "max" of all the k's tested so far
      }
      Mt=MaxType(Mn,longRunVar) #Setting the Max Type equal to Mt
    }
    if(Mt>critVal){ #Comparing the Max Type to Asymptotic value
      Cor=Cor+1 #Counting this as a change we detected
    }
    Mty=c(Mty, MaxType(Mn,longRunVar))
  }
  coverageProbability <- Cor/simDims[1]
  # Just returning the coverage probabilitiy for right now, may want to return other stuff later as well
  return(coverageProbability)
}

# Moved function declarations outside of loops
# Should probably talk to the group about this
TestStat<-function(y,x,n){ #Test Statistic (Zn) from Aue and ALexander
  (1/(sqrt(n)))*((sum(y[c(1,floor(n*x))]))-(((floor(n*x))/(n))*(sum(y[c(1,n)]))))
}

MaxType=function(m,w){ #Max Type Function from Aue Alexander
  (1/w)*m
}