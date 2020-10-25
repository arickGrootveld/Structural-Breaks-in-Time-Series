############################ indepDatagen.R ##################################
# Data generation function for Likelihood and CUSUM groups to test their algorithms
# against
#
# v2.0.0
# Contributors: Arick Grootveld
#############################################################################
indepDatagen <- function(simLen, numSims, param1=0.1, param2=0.2, randomBreaks=FALSE, breakLocations=5, distFamily="N", seed=0){
  ### Decided Simulation parameters ###
  
  # simLen <- 50                 # Length of each simulation (N)
  # 
  # numSims <- 1000               # Number of simulations in Monte Carlo process (M)
  # 
  # randomBreaks <- FALSE        # Boolean to decide if breaks locations will be 
  #                              # randomized
  # 
  # 
  # # If the random breaks is set to false, then 
  # # the location of the break will be set by this
  # # parameter. If this is a vector then it will
  # # cycle through break locations
  # breakLocations <- 0.5*simLen   
  #      
  # # Example values of breakLocations:
  # #         breakLocations <- c(3,5,2) # This will cause the break location to 
  # #                                    # cycle through these locations with 
  # #                                    # each new sequence of samples
  # #
  # #         breakLocations <- 5        # This will cause there to always be a 
  # #                                    # break at the 5th element of each 
  # #                                    # sequence
  # 
  # 
  # distFamily <- "N"            # String value from the set {"T", "N", "C", "E", "R"}, 
  #                              # corresponding to t, normal, chi-squared, and 
  #                              # exponential distributions, and R 
  #                              # corresponding to a random distribution 
  #                              # from this set. This value decides what family of 
  #                              # distribution will be used for the simulation
  # 
  # 
  # # parameters to modify between breaks
  # # Their exact meaning depends on the family of distribution, i.e.:
  # #   Normal:       parameters are the means of the before and after distributions
  # #   T:            parameters are the degrees of freedom of the two distributions (rounds up when non-integer value)
  # #   Chi-squared:  parameters are the degrees of freedom of the two distributions
  # #   Exponential:  parameters are the rate of the distributions (\lambda)
  # param1 <- 0.1
  # param2 <- 0.5
  # 
  # # TODO: allow user to specify that a random value should be selected from
  # # TODO: a range of values for the two parameters
  # 
  # # Seed used for random number generation (if set to >0, will use a random seed)
  # seed = 0
  
  ### Calculated System Parameters ###
  
  # Setting the seed of the simulation to make it recreatable
  if(seed > 0){
    set.seed(seed)
  }else if(exists('.Random.seed')){
    rm(.Random.seed, envir=globalenv())
  }
  
  # If dist family is "R", then it will select randomly from set of dist families
  if(distFamily == "R"){
    distFamily <- sample(c("N", "T", "C", "E"))
  }
  # TODO: make random dist family actually generate a 
  # TODO: vector of random distribution families so that each individual sim
  # TODO: in the Monte Carlo sim has a chance at any of the distributions
  
  # Round up on params for T distribution
  if(distFamily == "T"){
    param1 <- ceiling(param1)
    param2 <- ceiling(param2)
  }
  
  simMatrix <- matrix(0, numSims, simLen) # Preallocating the simulation matrix (MxN)
  
  # Conditional Structure to take care of assigning where breaks should happen
  if(randomBreaks){
    # Preallocating the break locations (cannot be either end of the simulation)
    # since that would result in a simulation with only one type of distribution
    breakLocs <- sample(2:(simLen - 1), numSims)
  } else {
    if(length(breakLocations) == 1){
      # If one break location was set then we simply say all sims should
      # have a break at that location
      breakLocs <- c(rep(breakLocations, numSims))
    }else{
      numBreakLocs <- length(breakLocations)
      if((numSims %% numBreakLocs) == 0){
        # If the break locations vector has a length that
        # evenly divides into the number of simulations then 
        # we simply repeat a proportional number of times
        numRepeats <- numSims / length(breakLocations)
        breakLocs <- c(rep(breakLocations, numRepeats))
      } else {
        # Other wise we repeat the break locations vector until we can't get
        # a full cycle in, and then fill it in with a partial vector
        numRepeats <- floor(numSims / length(breakLocations))
        partialIndex <- (numSims %% numBreakLocs)
        breakLocs <- c(rep(breakLocations, numRepeats), breakLocations[1:partialIndex])
      }
    }
  }
  
  # Assigning values based on the distribution family selected
  
  # Looping through all simulations to create rows of the simMatrix
  for (i in c(1:numSims)){
    if(distFamily == "N"){
      # Normal distributions
      simMatrix[i, 1:(breakLocs[i])] <- rnorm(length(1:(breakLocs[i])), mean=param1, sd=1)
      simMatrix[i, ((breakLocs[i]+1):simLen)] <- rnorm(length(((breakLocs[i]+1):simLen)), mean=param2, sd=1)
    }else if(distFamily == "T"){
      # T Distribution
      simMatrix[i, 1:(breakLocs[i])] <- rt(length(1:(breakLocs[i])), param1)
      simMatrix[i,((breakLocs[i]+1):simLen)] <- rt(length(((breakLocs[i]+1):simLen)), param2)
    }else if(distFamily == "C"){
      # Chi-squared dist
      simMatrix[i, 1:(breakLocs[i])] <- rchisq(length(1:(breakLocs[i])), param1)
      simMatrix[i, ((breakLocs[i]+1):simLen)] <- rchisq(length(((breakLocs[i]+1):simLen)), param2)
    }else if(distFamily == "E"){
      # Exponential distribution
      simMatrix[i, 1:(breakLocs[i])] <- rexp(length(1:(breakLocs[i])), rate=param1)
      simMatrix[i, ((breakLocs[i]+1):simLen)] <- rexp(length(((breakLocs[i]+1):simLen)), param2)
    }
  }
  return(simMatrix)
}
