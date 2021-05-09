############################ normalIndepMICBreakEstimator.R ############
# Function to calculate Modified Swartz Information Criterion, for detecting
# statistical breaks in normally distributed independent data
#
# v1.0.0
# Contributors: Andrea Scolari, Arick Grootveld, Ramadha Piyadi Gamage
#############################################################################
normIndepMICCalc <- function(simMatrix, alpha=0.05){
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
  # Critical value is calculated as the right-tailed value of a chi-squared dist.
  # based on theoretically derived values
  criticalValue <- qchisq(p=1-alpha, df=d)
  
  # For loop through all simulations
  for (simNum in c(1:simDims[1])){
    
    simSeq = simMatrix[simNum, ]
    
    # Filling the SIC1 array with INF values, so that we don't select any of the
    # default values when calculating the test statistic
    MIC <- rep(-Inf, length(simSeq))
    
    # The number of samples in the current simulation
    n = length(simSeq)
    
    
    # For loop through each index of potential break
    for (k in c(2:(n-2))){
      # Grabbing the sequence from before where we think the break happened
      before <- simSeq[1:k]
      # Grabbing the sequence from after the location the break happened
      after <- simSeq[(k+1):(n)]
      
      # MIC1[k] <- 2*n + 2*k * log(mean(before)) + 2*(n-i) * log(mean(after)) + (1+((2*i/n) - 1)^2) * log(n)
      MIC[k] <- 2* (n*log(var(simSeq)) - k*log(var(before)) - (n-k)*log(var(after))) - (((2*k/n) - 1)^2) * log(n) -(4 * log(n))
    }
    
    breakLoc <- which.max(MIC)
    
    testStatistic[1,simNum] <- max(MIC)
    
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