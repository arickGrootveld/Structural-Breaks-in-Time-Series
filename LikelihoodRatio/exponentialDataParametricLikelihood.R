############################ exponentialDataParametricLikelihoodCalc.R #######
# Function to calculate Likelihood Ratios, for detecting statistical breaks
# in exponentially distributed independent data
#
# v1.0.3
# Contributors: Andrea Scolari, Arick Grootveld, Ramadha Piyadi Gamage
#############################################################################
expIndepLRCalc <- function(simMatrix, alpha=0.05){
  # Calculated parameters to be used for the simulation
  
  simDims = dim(simMatrix)
  
  
  breakDetected = matrix(0,1, simDims[1])
  
  pVals <- matrix(0, 1, simDims[1])
  
  testStatistic <- matrix(0, 1, simDims[1])
  
  # Preallocating an array to contain the indexes of the estimated breaks, so we can calculate the break point location accuracy of
  # our algorithms
  detectedBreakIndexes = matrix(-1, 1, simDims[1])
  
  # Calculatings for the critical value, as provided by Ramadha Piyadi Gamage
  d <- 1 # d is the number of (unknown) parameters to be estimated
  a.n <- sqrt(2*log(log(simDims[2])))
  b.n <- a.n^2 + 0.5*log(log(log(simDims[2]))) - log(gamma(d/2))
  criticalValue <- ((-1/a.n)*log(log((1-alpha + exp(-2*exp(b.n)))^(-0.5))) + (b.n/a.n))^2 - log(simDims[2])
  
  # For loop through all simulations
  for (simNum in c(1:simDims[1])){
    
    logLikelihoodRatio <- double()
    simSeq = simMatrix[simNum, ]
    
    # The current simulation
    simSeqLen = length(simSeq)
    
    # denominators <- matrix(0, 1, length(simSeq))
    # 
    # # Calculating the numerator of the likelihood ratio as the pdf of each sample coming from a distribution with mean equal to the sample mean
    # numerator = prod(dexp(simSeq, rate=1/mean(simSeq)))
    
    # For loop through each index of potential break
    for (i in c(1:(simSeqLen-1))){
      # Grabbing the sequence from before where we think the break happened
      before <- simSeq[1:i]
      # Grabbing the sequence from after the location the break happened
      after <- simSeq[(i+1):(simSeqLen)]
      
      # Calculating the probablities of each of the samples being from the before distribution
      # beforePdf1 <- prod(dexp(before, rate=param1))
      # beforePdf2 <- prod(dexp(before, rate=param2))
      # 
      # # Calculating the probabilities of each of the samples from the after distribution
      # afterPdf1 <- prod(dexp(after, rate=param1))
      # afterPdf2 <- prod(dexp(after, rate=param2))
      # 
      # # Probability break was from dist1 to dist2, or from dist2 to dist1
      # option1 <- beforePdf1 * afterPdf2
      # option2 <- beforePdf2 * afterPdf1
      
      # Whichever one is larger of the posibilities is the denominator of the
      # likelihood ratio at this location
      #denominators[i] <- max(c(option1, option2))
      # denominators[i] <- option1
      
      logLikelihoodRatio[i] <- 2*((simSeqLen*log(mean(simSeq))) - (i*log(mean(before))) - (simSeqLen-i)*log(mean(after)))
    }
  
    breakLoc <- which.max(logLikelihoodRatio)
    
    testStatistic[1, simNum] <- max(logLikelihoodRatio)
    
    # If our test statistic is NaN, then we want to reject
    # This only happens if the maximal LR is negative
    if(is.na(testStatistic[1,simNum])){
      breakDetected[1,simNum] = 0
    }
    # If a break is detected, then we set the enum value to 1 (true)
    else if(testStatistic[1, simNum] > criticalValue){
      detectedBreakIndexes[1,simNum] = breakLoc
      breakDetected[1,simNum] = 1
    }
    # Technically redundant, since its already defaulted to zeros, but just to be safe
    else{
      breakDetected[1,simNum] = 0
    }
  }
  
  # Calculating the coverage probability
  numDetected <- sum(breakDetected)
  
  coverageProbability <- numDetected / simDims[1]
  return(c(coverageProbability, detectedBreakIndexes))
}