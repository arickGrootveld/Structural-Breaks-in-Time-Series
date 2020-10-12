# Parameters to be used for the simulation
param1 <- 0.1
param2 <- 1

alpha = 0.05

simDims = dim(simMatrix)
denominators <- matrix(0, simDims)

breakDetected = matrix(0,1, simDims[1])

pVals <- matrix(0, 1, simDims[1])

for (simNum in c(1:simDims[1])){
  simSeq = simMatrix[simNum, ]
  
  denominators <- matrix(0, 1, length(simSeq))
  
  # Calculating pdfs of each sample coming from each parameter only
  fromDist1Prob = prod(dnorm(simSeq, mean=param1, sd=1))
  fromDist2Prob = prod(dnorm(simSeq, mean=param2, sd=1))
  
  numerator = max(c(fromDist1Prob, fromDist2Prob))
  
  for (i in c(2:length(simSeq)-1)){
    # Grabbing the sequence from before where we think the break happened
    before <- simSeq[1:i]
    # Grabbing the sequence from after the location the break happened
    after <- simSeq[i:length(simSeq)]
    
    # Calculating the probablities of each of the samples being from the before distribution
    beforePdf1 <- prod(dnorm(before, mean=param1, sd=1))
    beforePdf2 <- prod(dnorm(before, mean=param2, sd=1))
    # Calculating the probabilities of each of the samples from the after distribution
    afterPdf1 <- prod(dnorm(after, mean=param1, sd=1))
    afterPdf2 <- prod(dnorm(after, mean=param2, sd=1))
    
    # Probability break was from dist1 to dist2, or from dist2 to dist1
    option1 <- beforePdf1 * afterPdf2
    option2 <- beforePdf2 * afterPdf1
    
    # Whichever one is larger of the posibilities is the likelihood at this
    # location
    denominators[i] <- max(c(option1, option2))
  }
  
  breakLoc <- which.max(denominators[1:length(denominators)-1])
  denominatorValue <- max(denominators[1:length(denominators)-1])
  
  likelihoodRatio = -2 * log(numerator / denominatorValue)
  
  pVals[1,simNum] = 1 - pchisq(likelihoodRatio, df=1)
  # If a break is detected, then we set the enum value to 1 (true)
  if(pVals[1, simNum] < alpha){
    breakDetected[1,simNum] = 1
  }
}

# Calculating the coverage probability
numDetected <- sum(breakDetected)

coverageProbability <- numDetected / simDims[1]