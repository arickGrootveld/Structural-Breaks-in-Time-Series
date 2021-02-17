############################ CUSUMCalc.R ##################################
# Function to calculate Likelihood Ratios, for detecting statistical breaks
# in normally distributed independent data
#
# v1.0.2
# Contributors: Andrea Scolari, Arick Grootveld
#############################################################################
normalIndepLRCalc <- function(param1, param2, simMatrix, alpha=0.05){
  # Calculated parameters to be used for the simulation
  criticalValue = -log(-0.5*log(1-alpha))
  
  simDims = dim(simMatrix)
  denominators <- matrix(0, simDims+1)
  
  breakDetected = matrix(0,1, simDims[1])
  
  pVals <- matrix(0, 1, simDims[1])
  
  testStatistic <- matrix(0, 1, simDims[1])
  
  for (simNum in c(1:simDims[1])){
    simSeq = simMatrix[simNum, ]
    
    simSeqLen = length(simSeq)
    
    denominators <- matrix(0, 1, length(simSeq))
    
    # Calculating the numerator of the likelihood ratio as the pdf of each sample coming from a distribution with mean equal to the sample mean
    numerator = prod(dnorm(simSeq, mean=mean(simSeq), sd=1))
    
    for (i in c(1:(simSeqLen-1))){
      # Grabbing the sequence from before where we think the break happened
      before <- simSeq[1:i]
      # Grabbing the sequence from after the location the break happened
      after <- simSeq[(i+1):(simSeqLen)]
      
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
    # denominators[i+1] = numerator
    
    breakLoc <- which.max(denominators)
    
    #### End of test lines
    likelihoodRatio = max(-2 * log(numerator / denominators))
    
    # Special value that Ramadha is having us calculate
    testStatistic[1,simNum] = (2*log(log(simDims[2])))^(1/2) * (likelihoodRatio)^(1/2) - (2*log(log(simDims[2])) + log(log(log(simDims[2]))))
    
    pVals[1,simNum] = 1 - pchisq(likelihoodRatio, df=1)
    
    #if(pVals[1, simNum] < alpha){
    #  breakDetected[1, simNum] = 1
    #}
    
    # If our test statistic is NaN, then we want to reject
    # This only happens if the maximal LR is negative
    if(is.na(testStatistic[1,simNum])){
      breakDetected[1,simNum] = 0
    }
    # If a break is detected, then we set the enum value to 1 (true)
    else if(testStatistic[1, simNum] > criticalValue){
      breakDetected[1,simNum] = 1
    }
    # Technically redundant, since its already defaulted to zeros, but just to be safe
    else{
      breakDetected[1,simNum] = 0
    }
  }
  
  # Calculating the coverage probability
  numDetected <- sum(breakDetected)
  
  return(coverageProbability <- numDetected / simDims[1])
}