############################ timeSeriesDatagen.R ##################################
# Script to generate AR data for Time Series Change Point Simulations
#
# v1.0.0
# Contributors: Arick Grootveld
#############################################################################

timeSeriesDatagen <- function(simLen, numSims, param1, param2, randomBreaks = FALSE, breakLocations = 5, processName = 'AR', seed=0, noiseDist = "N", noiseVar = 1, startingPosition = 0, randomStart=TRUE){

  # Setting the seed of the simulation to make it recreatable
  if(seed > 0){
    set.seed(seed)
  }else if(exists('.Random.seed')){
    rm(.Random.seed, envir=globalenv())
  }

  # Simulation code written for AR Process, eventually will add support for GARCH, ARMA, etc.
  if(processName == "AR"){
    if(breakLocations > simLen){
      badParamsErrorString <- paste('The break location (', toString(breakLocations), ') is larger than the total simulation length (', toString(simLen), ')')
      stop(simpleError(badParamsErrorString))
    }

    ## TODO: Implement these at some point
    if(randomBreaks==TRUE){
      notImplementedErrorString <- "Random breaks not implemented yet"
      stop(simpleError(notImplementedErrorString))
    }
    if(length(breakLocations) > 1){
      notImplementedErrorString <- "More than one potential break location not implemented yet"
      stop(simpleError(notImplementedErrorString))
    }
    if(noiseDist != "N"){
      notImplementedErrorString <- "Non-normal noise terms have not been implemented yet"
      stop(simpleError(notImplementedErrorString))
    }

    ## If all params pass basic sanity checks we go ahead and implement the AR Process
    # Generating the noise terms that will be added to the AR process values at each step
    noiseTerms <- rnorm(n=simLen*numSims, mean=0, sd = sqrt(noiseVar))

    simMatrix <- matrix(NA, numSims, simLen) # Preallocating the simulation matrix (MxN)

    # The first sample needs to have a variance of \sigma_y^2 
    # based on the assumptions made on page 173 (pdf version) 
    # of Art B. Owens "Empirical Likelihood" book. This is 
    # the assumption described in the text between equations 
    # 8.2 and 8.3 
    if(randomStart == TRUE){
      # simMatrix[,1] <- startingPosition + c(noiseTerms[seq(from=1, to=length(noiseTerms), by=simLen)])
      simMatrix[,1] <- startingPosition + rnorm(n=numSims, mean=0, sd = sqrt(noiseVar/(1 - param1^2)))
    }else{
      simMatrix[,1] <- startingPosition
    }
    # simMatrix[,1] <- c(noiseTerms[seq(from=1, to=length(noiseTerms), by=simLen)])


    # Looping through each simulation and calculating each of the values for each of them
    for(i in 1:numSims){
      # Looping through each element of each simulation and calculating the value from the noise
      # and previous terms
      for(j in 2:simLen){
        # If the break hasn't happened yet we use the first param
        if(j < breakLocations){
          simMatrix[i, j] <- param1 * simMatrix[i,(j-1)] + noiseTerms[j + (simLen*(i-1))]
        }else{
          simMatrix[i,j] <- param2 * simMatrix[i,(j-1)] + noiseTerms[j + (simLen*(i-1))]
        }

      }
    }

    return(simMatrix)

  # Otherwise we will throw an error because the specified process name is not on
  # the list of known processes
  }else{
    badParamsErrorString <- paste('The parameter: [', processName, "] is not a known process name. Known names include 'AR'")
    badPassedParamError <- simpleError(badParamsErrorString)
    stop(badPassedParamError)
  }

}