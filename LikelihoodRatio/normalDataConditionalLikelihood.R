############################ normalDataConditionalLikelihood.R ##################################
# Function to calculate Conditional Likelihood Ratio, for detecting statistical breaks
# in normally distributed time series data suspected to come from a AR process
#
# The conditional likelihood implementation is taken from 
# "Empirical Likelihood" by Art B. Owen. Specifically from the equation
# on page 173 (pdf version) equation 8.5
# 
#
# v1.0.0
# Contributors: Andrea Scolari, Arick Grootveld, Ramadha Piyadi Gamage
#
# Inputs: 
#         simMatrix: matrix of data to use for simulations (generally as output 
#                    of indepDataGen script). Each row should be a sequence of 
#                    independently distributed data, such that the first column
#                    contains the first sample of each sequence
#         mu:     The mean of the noise term for the AR process
#         sigma:  The standard deviation of the noise term for the AR process
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
normalCondLRCalc <- function(simMatrix, mu, sigma, alpha=0.05){
  # Calculated parameters to be used for the simulation
  criticalValue = -log(-0.5*log(1-alpha))
  
  simDims = dim(simMatrix)
  denominators <- matrix(0, simDims+1)
  
  breakDetected = matrix(0,1, simDims[1])
  
  pVals <- matrix(0, 1, simDims[1])
  
  testStatistic <- matrix(0, 1, simDims[1])
  
  # Preallocating an array to contain the indexes of the estimated breaks, so we can calculate the break point location accuracy of
  # our algorithms
  detectedBreakIndexes = matrix(-1, 1, simDims[1])
  
  # For loop through all simulations
  for (simNum in c(1:simDims[1])){
    simSeq = simMatrix[simNum, ]
    
    # The current simulation
    simSeqLen = length(simSeq)
    
    denominators <- matrix(0, 1, length(simSeq))
    
    # Estimating the AR process parameter using the 
    # OLS estimator as shown in this video:
    # https://www.youtube.com/watch?v=YvB3k00SF7w
    # OLS_est_of_AR_param <- (sum(simSeq[2:simSeqLen] * simSeq[1:(simSeqLen - 1)])) / (sum(simSeq[1:(simSeqLen-1)]^2))
    
    # estOfARParam <- as.numeric(arima(simSeq, order=c(1,0,0), method="ML")$coef[1])
    estOfARParam <- as.numeric(ar(simSeq, order.max=1, method='yw', aic=FALSE)$ar)
    
    # Estimating the numerator portion from the conditional likelihood
    numerator <- condLikelihood(timeSeries = simSeq, mu = mu, sigma = sigma, beta_1 = estOfARParam)
    
    # For loop through each index of potential break
    for (i in c(3:(simSeqLen-3))){
      # Grabbing the sequence from before where we think the break happened
      before <- simSeq[1:i]
      # Grabbing the sequence from after the location the break happened
      after <- simSeq[(i+1):(simSeqLen)]
      
      # beforeEstOfARParam <- as.numeric(arima(before, order=c(1,0,0), method="ML")$coef[1])
      # afterEstOfARParam <- as.numeric(arima(after, order=c(1,0,0), method="ML")$coef[1])
      beforeEstOfARParam <-  as.numeric(ar(before, order.max=1, method='yw', aic=FALSE)$ar)
      afterEstOfARParam <-  as.numeric(ar(after, order.max=1, method='yw', aic=FALSE)$ar)
      
      # Calculating the probabilities of each of the samples being from the before distribution
      beforePdf <- condLikelihood(timeSeries = before, mu=mu, sigma=sigma, beta_1 = beforeEstOfARParam)
      
      # Calculating the probabilities of each of the samples from the after distribution
      afterPdf <- condLikelihood(timeSeries = after, mu=mu, sigma = sigma, beta_1 = afterEstOfARParam)
      
      denominators[i] <- beforePdf * afterPdf
    }
    breakLoc <- which.max(denominators)
    # denominators[which(denominators == 0)] = Inf
    
    likelihoodRatio = max(-2 * log(numerator / denominators))
    
    # Special value that Ramadha is having us calculate
    testStatistic[1,simNum] = (2*log(log(simDims[2])))^(1/2) * (likelihoodRatio)^(1/2) - (2*log(log(simDims[2])) + log(log(log(simDims[2]))))
    
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
  
  detectedBreakProportion <- numDetected / simDims[1]
  # Returning the proportion of detected breaks to total sample size as the 
  # first returned value, and the detected break locations as an array after 
  # the first value. Any values with -1 mean that no break was detected for that
  # simulation
  return(c(detectedBreakProportion, detectedBreakIndexes))
}


# Function for calculating the conditional likelihood of a series of samples, assuming the samples
# come from an AR(1) process with a noise term ~ N(mu, sigma^2), and with AR param beta_1
condLikelihood <- function(timeSeries, mu, sigma, beta_1){
  numEls <- length(timeSeries)
  likelihoodHolder <- 1
  
  for(i in 2:numEls){
    likelihoodHolder <- likelihoodHolder * (1/(sqrt(2*pi) * sigma)) * exp(-(1/(2*sigma^2)) * ((timeSeries[i] - mu) - (beta_1 * (timeSeries[i - 1] - mu)) )^2 )
  }
  return(likelihoodHolder)
}