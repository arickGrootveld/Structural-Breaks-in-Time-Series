############################ exponentialDataParametricLikelihoodCalc.R #######
# Function to calculate Likelihood Ratios, for detecting statistical breaks
# in exponentially distributed independent data
#
# v1.0.3
# Contributors: Andrea Scolari, Arick Grootveld, Ramadha Piyadi Gamage
#
# Inputs: 
#         simMatrix: matrix of data to use for simulations (generally as output 
#                    of indepDataGen script). Each row should be a sequence of 
#                    independently distributed data, such that the first column
#                    contains the first sample of each sequence
#         alpha: the significance level for the estimator
# Outputs:
#         outputVals: array of values containing all the results from the 
#                     simulation
#             - outputVals[1]: the proportion of the detected breaks from
#                              the number of sequences that were input
#             - outputVals[2:-1]: the indexes of the detected breaks. If a value
#                                 is -1, then no break was detected for that 
#                                 sequence
#
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
    
    # For loop through each index of potential break
    for (i in c(1:(simSeqLen-1))){
      # Grabbing the sequence from before where we think the break happened
      before <- simSeq[1:i]
      # Grabbing the sequence from after the location the break happened
      after <- simSeq[(i+1):(simSeqLen)]
      
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
  
  # Calculating the proportion of detected breaks
  numDetected <- sum(breakDetected)
  detectedBreakProportion <- numDetected / simDims[1]
  # Returning the proportion of detected breaks to total sample size as the 
  # first returned value, and the detected break locations as an array after 
  # the first value. Any values with -1 mean that no break was detected for that
  # simulation
  return(c(detectedBreakProportion, detectedBreakIndexes))
}