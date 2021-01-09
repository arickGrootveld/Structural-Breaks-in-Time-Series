# Parameters to be used for the simulation
param1 <- 0.1
param2 <- 0.5

alpha = 0.05
criticalValue = -log(-0.5*log(1-alpha))

simDims = dim(simMatrix)
denominators <- matrix(0, simDims+1)

breakDetected = matrix(0,1, simDims[1])

pVals <- matrix(0, 1, simDims[1])

testStatistic <- matrix(0, 1, simDims[1])

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
    
    # Whichever one is larger of the posibilities is the denominator of the
    # likelihood ratio at this location
    denominators[i] <- max(c(option1, option2))
  }
  
  # Very important step
  # We add the possibility of no breaks into the denominator, so the numberator becomes a subset of the denominator
  # This is part of the requirement for Wilks's theorm to be true
  denominators[i+1] = numerator
  
  breakLoc <- which.max(denominators)
  
  #### The Values contained between these lines are for test purposes only, if you see this, please delete them
  testValues = -2 * log(numerator / denominators)
  
  #### End of test lines
  likelihoodRatio = max(-2 * log(numerator / denominators))
  
  # Special value that Ramadha is having us calculate
  testStatistic[1,simNum] = (2*log(log(simDims[2])))^(1/2) * (likelihoodRatio)^(1/2) - (2*log(log(simDims[2])) + log(log(log(simDims[2]))))
  
  pVals[1,simNum] = 1 - pchisq(likelihoodRatio, df=1)
  # If a break is detected, then we set the enum value to 1 (true)
  if(testStatistic[1, simNum] > criticalValue){
    breakDetected[1,simNum] = 1
  }
  
  #if(pVals[1, simNum] < alpha){
  #  breakDetected[1,simNum] = 1
  #}
}

# Calculating the coverage probability
numDetected <- sum(breakDetected)

coverageProbability <- numDetected / simDims[1]