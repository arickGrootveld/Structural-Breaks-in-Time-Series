
## This script generates multiple datasets of time series (AR(1))
## data, with breaks at specified locations. The script
## computes the ability of three methods to detect the presence of these
## breaks: Likelihood Ratio, and SIC.

## This code is effectively very similar to the runNormalDataProcessAndGatherResults
## code, but with the exception that we use exponentially distributed data 
## instead of normally distributed data for generating the samples

# Importing datagen function
source('indepDatagen.R')

# Importing the Likelihood Ration function
source('LikelihoodRatio/exponentialDataParametricLikelihood.R')

# Importing SIC and MIC calculation functions
source('SIC/exponentialIndepSICBreakEstimator.R')
source('SIC/exponentialIndepMICBreakEstimator.R')

# Initializing variables

# Modifiable Parameters
################################################################################
# Difference in mean between distributions
param1 = 1
param2 = 2
# Number of simulations to perform
numSims = 100

# Seed for simulations (if set to 0, the seed is random, and any previously 
# set seed is globally cleared)
curSeed = 0

## Parameters for the Likelihood ratio, SIC and MIC algs.
significanceLevel = 0.05
# Whether to run LR calculation in this expirement or not
runLR = 1

# Whether to run SIC calculations in this experiment or not
runSIC = 1

# Whether to run MIC calculations in this experiment or not
runMIC = 1

## Parameters defining the type of output the user wants to see
# This will tell the code to either output the percent error in detected break locations, or not (1 to run, 0 to not)
showBreakLocDetectAcc = 0

# This will tell the code to output coverage probability curves for the method (TODO: Implement this)
showCoverageProbabilityCurves = 0
# If the above variable is 1, then this variable will set the step size of the variable used to linearly interpolate the 
# line of the coverage probability curve(as a fraction percent out of 1, i.e.
# 0.1 would plot the score in intervals of 10% of the total sequence length)
covProbParamStepSize = 0.04

# This variables dictates whether a parameter sweep graph is generated as 
# per page 163 of the paper: "The Likelihood Ratio Test for the change point 
# problem with exponentially distributed random variables"
showParamSweepPlot = 0
# Parameter to determine if the parameter sweep plot should be saved or not
saveParamSweepPlot = 0
# Parameter to determine the number of steps of the parameter sweep plot
paramSweepPlotStepNumber = 100

# Variable declarations
################################################################################
## Arrays to store results of likelihood ratio and SIC calculations
nEquals30LR = c(0,0,0,0)
nEquals100LR = c(0,0,0,0)
nEquals200LR = c(0,0,0,0)

nEquals30SIC = c(0,0,0,0)
nEquals100SIC = c(0,0,0,0)
nEquals200SIC = c(0,0,0,0)

nEquals30MIC = c(0,0,0,0)
nEquals100MIC = c(0,0,0,0)
nEquals200MIC = c(0,0,0,0)

## Variables to store the average break location detected for each technique
breakIndexes30LR = c(0, 0, 0, 0)
breakIndexes100LR = c(0, 0, 0, 0)
breakIndexes200LR = c(0, 0, 0, 0) 

breakIndexes30SIC = c(0,0,0,0)
breakIndexes100SIC = c(0,0,0,0)
breakIndexes200SIC = c(0,0,0,0)

breakIndexes30MIC = c(0,0,0,0)
breakIndexes100MIC = c(0,0,0,0)
breakIndexes200MIC = c(0,0,0,0)


# Performing simulations
################################################################################
## n = 30
simLen <- 30
### K=0.2n
breakLocs <- 0.2 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed, distFamily='E')
# LR Calculation
if(runLR==1){
  retVal = expIndepLRCalc(simMatrix, alpha=significanceLevel)
  nEquals30LR[1] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes30LR[1] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes30LR[1] = -1
    }
  }
}

if(runSIC==1){
  # SIC Calc
  retVal = expIndepSICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals30SIC[1] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes30SIC[1] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes30SIC[1] = -1
    }
  }
}

if(runMIC==1){
  # MIC Calc
  retVal = expIndepMICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals30MIC[1] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes30MIC[1] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes30MIC[1] = -1
    }
  }
}

### K=0.3n
breakLocs <- 0.3 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed, distFamily='E')
# LR Calculation

if(runLR==1){
  retVal = expIndepLRCalc(simMatrix, alpha=significanceLevel)
  nEquals30LR[2] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes30LR[2] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes30LR[2] = -1
    }
  }
}

if(runSIC==1){
  # SIC Calc
  retVal = expIndepSICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals30SIC[2] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes30SIC[2] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes30SIC[2] = -1
    }
  }
}


if(runMIC==1){
  # SIC Calc
  retVal = expIndepMICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals30MIC[2] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes30MIC[2] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes30MIC[2] = -1
    }
  }
}

### K=0.5n
breakLocs <- 0.5 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed, distFamily='E')

# LR Calculation
if(runLR==1){
  retVal = expIndepLRCalc(simMatrix, alpha=significanceLevel)
  nEquals30LR[3] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes30LR[3] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes30LR[3] = -1
    }
  }
}

if(runSIC==1){
  # SIC Calc
  retVal = expIndepSICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals30SIC[3] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes30SIC[3] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes30SIC[3] = -1
    }
  }
}


if(runMIC==1){
  # MIC Calc
  retVal = expIndepMICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals30MIC[3] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes30MIC[3] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes30MIC[3] = -1
    }
  }
}

### K=0.75n
breakLocs <- 0.75 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed, distFamily='E')
# LR Calculation
if(runLR==1){
  retVal = expIndepLRCalc(simMatrix, alpha=significanceLevel)
  nEquals30LR[4] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes30LR[4] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes30LR[4] = -1
    }
  }
}

if(runSIC==1){
  # SIC Calc
  retVal = expIndepSICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals30SIC[4] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Getting all break detection indexes so we can calculate the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes30SIC[4] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes30SIC[4] = -1
    }
  }
}


if(runMIC==1){
  # MIC Calc
  retVal = expIndepMICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals30MIC[4] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes30MIC[4] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes30MIC[4] = -1
    }
  }
}

################################################################################
## n = 100
simLen <- 100
### K=0.2n
breakLocs <- 0.2 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed, distFamily='E')
# LR Calculation

if(runLR==1){
  retVal = expIndepLRCalc(simMatrix, alpha=significanceLevel)
  nEquals100LR[1] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes100LR[1] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes100LR[1] = -1
    }
  }
}

if(runSIC==1){
  # SIC Calc
  retVal = expIndepSICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals100SIC[1] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes100SIC[1] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes100SIC[1] = -1
    }
  }
}

if(runMIC==1){
  # MIC Calc
  retVal = expIndepMICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals100MIC[1] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes100MIC[1] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes100MIC[1] = -1
    }
  }
}

### K=0.3n
breakLocs <- 0.3 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed, distFamily='E')
# LR Calculation

if(runLR==1){
  retVal = expIndepLRCalc(simMatrix, alpha=significanceLevel)
  nEquals100LR[2] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes100LR[2] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes100LR[2] = -1
    }
  }
}

if(runSIC==1){
  # SIC Calc
  retVal = expIndepSICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals100SIC[2] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes100SIC[2] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes100SIC[2] = -1
    }
  }
}


if(runMIC==1){
  # MIC Calc
  retVal = expIndepMICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals100MIC[2] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes100MIC[2] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes100MIC[2] = -1
    }
  }
}

### K=0.5n
breakLocs <- 0.5 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed, distFamily='E')
# LR Calculation

if(runLR==1){
  retVal = expIndepLRCalc(simMatrix, alpha=significanceLevel)
  nEquals100LR[3] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes100LR[3] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes100LR[3] = -1
    }
  }
}

if(runSIC==1){
  # SIC Calc
  retVal = expIndepSICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals100SIC[3] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes100SIC[3] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes100SIC[3] = -1
    }
  }
}


if(runMIC==1){
  # MIC Calc
  retVal = expIndepMICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals100MIC[3] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes100MIC[3] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes100MIC[3] = -1
    }
  }
}

### K=0.75n
breakLocs <- 0.75 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed, distFamily='E')
# LR Calculation

if(runLR==1){
  retVal = expIndepLRCalc(simMatrix, alpha=significanceLevel)
  nEquals100LR[4] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes100LR[4] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes100LR[4] = -1
    }
  }
}

if(runSIC==1){
  # SIC Calc
  retVal = expIndepSICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals100SIC[4] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes100SIC[4] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes100SIC[4] = -1
    }
  }
}


if(runMIC==1){
  # MIC Calc
  retVal = expIndepMICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals100MIC[4] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes100MIC[4] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes100MIC[4] = -1
    }
  }
}

################################################################################
## n = 200
simLen <- 200
### K=0.2n
breakLocs <- 0.2 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed, distFamily='E')
# LR Calculation

if(runLR==1){
  retVal = expIndepLRCalc(simMatrix, alpha=significanceLevel)
  nEquals200LR[1] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes200LR[1] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes200LR[1] = -1
    }
  }
}

if(runSIC==1){
  # SIC Calc
  retVal = expIndepSICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals200SIC[1] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes200SIC[1] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes200SIC[1] = -1
    }
  }
}


if(runMIC==1){
  # MIC Calc
  retVal = expIndepMICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals200MIC[1] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes200MIC[1] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes200MIC[1] = -1
    }
  }
}


### K=0.3n
breakLocs <- 0.3 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed, distFamily='E')
# LR Calculation

if(runLR==1){
  retVal = expIndepLRCalc(simMatrix, alpha=significanceLevel)
  nEquals200LR[2] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes200LR[2] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes200LR[2] = -1
    }
  }
}

if(runSIC==1){
  # SIC Calc
  retVal = expIndepSICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals200SIC[2] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes200SIC[2] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes200SIC[2] = -1
    }
  }
}


if(runMIC==1){
  # MIC Calc
  retVal = expIndepMICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals200MIC[2] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes200MIC[2] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes200MIC[2] = -1
    }
  }
}

### K=0.5n
breakLocs <- 0.5 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed, distFamily='E')
# LR Calculation

if(runLR==1){
  retVal = expIndepLRCalc(simMatrix, alpha=significanceLevel)
  nEquals200LR[3] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes200LR[3] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes200LR[3] = -1
    }
  }
}


if(runSIC==1){
  # SIC Calc
  retVal = expIndepSICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals200SIC[3] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes200SIC[3] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes200SIC[3] = -1
    }
  }
}


if(runMIC==1){
  # MIC Calc
  retVal = expIndepMICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals200MIC[3] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes200MIC[3] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes200MIC[3] = -1
    }
  }
}

### K=0.75n
breakLocs <- 0.75 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed, distFamily='E')
# LR Calculation

if(runLR==1){
  retVal = expIndepLRCalc(simMatrix, alpha=significanceLevel)
  nEquals200LR[4] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes200LR[4] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes200LR[4] = -1
    }
  }
}


if(runSIC==1){
  # SIC Calc
  retVal = expIndepSICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals200SIC[4] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes200SIC[4] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes200SIC[4] = -1
    }
  }
}


if(runMIC==1){
  # MIC Calc
  retVal = expIndepMICCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals200MIC[4] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes200MIC[4] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes200MIC[4] = -1
    }
  }
}

################################################################################
## Constructing the output table for Likelihood Ratio and SIC
results.dataLR <- data.frame(
  breakLocs = c("K=0.2n", "K=0.3n", "K=0.5n", "K=0.75n"),
  nEquals30 = nEquals30LR,
  nEquals100 = nEquals100LR,
  nEquals200 = nEquals200LR
)

results.dataSIC <- data.frame(
  breakLocs = c("K=0.2n", "K=0.3n", "K=0.5n", "K=0.75n"),
  nEquals30 = nEquals30SIC,
  nEquals100 = nEquals100SIC,
  nEquals200 = nEquals200SIC
)


results.dataMIC <- data.frame(
  breakLocs = c("K=0.2n", "K=0.3n", "K=0.5n", "K=0.75n"),
  nEquals30 = nEquals30MIC,
  nEquals100 = nEquals100MIC,
  nEquals200 = nEquals200MIC
)

## Constructing the output table for the detected break locations
# Getting actual break locations so we can calculate error of the detected break locations
if(showBreakLocDetectAcc == 1){
  bAt30For_0_2n = round(30 * 0.2)
  bAt30For_0_3n = round(30 * 0.3)
  bAt30For_0_5n = round(30 * 0.5)
  bAt30For_0_75n = round(30 * 0.75)
  
  bAt100For_0_2n = round(100 * 0.2)
  bAt100For_0_3n = round(100 * 0.3)
  bAt100For_0_5n = round(100 * 0.5)
  bAt100For_0_75n = round(100 * 0.75)
  
  bAt200For_0_2n = round(200 * 0.2)
  bAt200For_0_3n = round(200 * 0.3)
  bAt200For_0_5n = round(200 * 0.5)
  bAt200For_0_75n = round(200 * 0.75)
  
  
  ## Calculating the error of the LR method's break location detection
  breakResults30_LR = c(-1,-1,-1,-1)
  breakResults100_LR = c(-1, -1, -1, -1)
  breakResults200_LR = c(-1,-1,-1,-1)
  # Calculating the error when n=30 (in terms of percent of total samples)
  # If a break was detected we can do the calculation, otherwise we leave it as -1
  if(breakIndexes30LR[1] > 0){
    breakResults30_LR[1] = 100 * abs(breakIndexes30LR[1] - bAt30For_0_2n) / 30
  }
  
  if(breakIndexes30LR[2] > 0){
    breakResults30_LR[2] = 100 * abs(breakIndexes30LR[1] - bAt30For_0_3n) / 30
  }
  
  if(breakIndexes30LR[3] > 0){
    breakResults30_LR[3] = 100 * abs(breakIndexes30LR[3] - bAt30For_0_5n) / 30
  }
  
  if(breakIndexes30LR[4] > 0){
    breakResults30_LR[4] = 100 * abs(breakIndexes30LR[4] - bAt30For_0_75n) / 30
  }
  
  # Calculating the error when n=100 (in terms of percent of total samples)
  if(breakIndexes100LR[1] > 0){
    breakResults100_LR[1] = 100 * abs(breakIndexes100LR[1] - bAt100For_0_2n) / 100
  }
  
  if(breakIndexes100LR[2] > 0){
    breakResults100_LR[2] = 100 * abs(breakIndexes100LR[2] - bAt100For_0_3n) / 100
  }
  
  if(breakIndexes100LR[3] > 0){
    breakResults100_LR[3] = 100 * abs(breakIndexes100LR[3] - bAt100For_0_5n) / 100
  }
  
  if(breakIndexes100LR[4] > 0){
    breakResults100_LR[4] = 100 * abs(breakIndexes100LR[4] - bAt100For_0_75n) / 100
  }
  
  # Calculating the error when n=200 (in terms of percent of total samples)
  if(breakIndexes200LR[1] > 0){
    breakResults200_LR[1] = 100 * abs(breakIndexes200LR[1] - bAt200For_0_2n) / 200
  }
  
  if(breakIndexes200LR[2] > 0){
    breakResults200_LR[2] = 100 * abs(breakIndexes200LR[2] - bAt200For_0_3n) / 200
  }
  
  if(breakIndexes200LR[3] > 0){
    breakResults200_LR[3] = 100 * abs(breakIndexes200LR[3] - bAt200For_0_5n) / 200
  }
  
  if(breakIndexes200LR[4] > 0){
    breakResults200_LR[4] = 100 * abs(breakIndexes200LR[4] - bAt200For_0_75n) / 200
  }
  
  breakResults.dataLR <- data.frame(
    breakLocs = c("K=0.2n", "K=0.3n", "K=0.5n", "K=0.75n"),
    nEquals30 = breakResults30_LR,
    nEquals100 = breakResults100_LR,
    nEquals200 = breakResults200_LR
  )
  
  
  ## Calculating the error of the SIC's break location detection
  breakResults30_SIC = c(-1,-1,-1,-1)
  breakResults100_SIC = c(-1, -1, -1, -1)
  breakResults200_SIC = c(-1,-1,-1,-1)
  # Calculating the error when n=30 (in terms of percent of total samples)
  # If a break was detected we can do the calculation, otherwise we leave it as -1
  if(breakIndexes30SIC[1] > 0){
    breakResults30_SIC[1] = 100 * abs(breakIndexes30SIC[1] - bAt30For_0_2n) / 30
  }
  
  if(breakIndexes30SIC[2] > 0){
    breakResults30_SIC[2] = 100 * abs(breakIndexes30SIC[1] - bAt30For_0_3n) / 30
  }
  
  if(breakIndexes30SIC[3] > 0){
    breakResults30_SIC[3] = 100 * abs(breakIndexes30SIC[3] - bAt30For_0_5n) / 30
  }
  
  if(breakIndexes30SIC[4] > 0){
    breakResults30_SIC[4] = 100 * abs(breakIndexes30SIC[4] - bAt30For_0_75n) / 30
  }
  
  # Calculating the error when n=100 (in terms of percent of total samples)
  if(breakIndexes100SIC[1] > 0){
    breakResults100_SIC[1] = 100 * abs(breakIndexes100SIC[1] - bAt100For_0_2n) / 100
  }
  
  if(breakIndexes100SIC[2] > 0){
    breakResults100_SIC[2] = 100 * abs(breakIndexes100SIC[2] - bAt100For_0_3n) / 100
  }
  
  if(breakIndexes100SIC[3] > 0){
    breakResults100_SIC[3] = 100 * abs(breakIndexes100SIC[3] - bAt100For_0_5n) / 100
  }
  
  if(breakIndexes100SIC[4] > 0){
    breakResults100_SIC[4] = 100 * abs(breakIndexes100SIC[4] - bAt100For_0_75n) / 100
  }
  
  # Calculating the error when n=200 (in terms of percent of total samples)
  if(breakIndexes200SIC[1] > 0){
    breakResults200_SIC[1] = 100 * abs(breakIndexes200SIC[1] - bAt200For_0_2n) / 200
  }
  
  if(breakIndexes200SIC[2] > 0){
    breakResults200_SIC[2] = 100 * abs(breakIndexes200SIC[2] - bAt200For_0_3n) / 200
  }
  
  if(breakIndexes200SIC[3] > 0){
    breakResults200_SIC[3] = 100 * abs(breakIndexes200SIC[3] - bAt200For_0_5n) / 200
  }
  
  if(breakIndexes200SIC[4] > 0){
    breakResults200_SIC[4] = 100 * abs(breakIndexes200SIC[4] - bAt200For_0_75n) / 200
  }
  
  breakResults.dataSIC <- data.frame(
    breakLocs = c("K=0.2n", "K=0.3n", "K=0.5n", "K=0.75n"),
    nEquals30 = breakResults30_SIC,
    nEquals100 = breakResults100_SIC,
    nEquals200 = breakResults200_SIC
  )
  
  
  ## Calculating the error of the MIC's break location detection
  breakResults30_MIC = c(-1,-1,-1,-1)
  breakResults100_MIC = c(-1, -1, -1, -1)
  breakResults200_MIC = c(-1,-1,-1,-1)
  # Calculating the error when n=30 (in terms of percent of total samples)
  # If a break was detected we can do the calculation, otherwise we leave it as -1
  if(breakIndexes30MIC[1] > 0){
    breakResults30_MIC[1] = 100 * abs(breakIndexes30MIC[1] - bAt30For_0_2n) / 30
  }
  
  if(breakIndexes30MIC[2] > 0){
    breakResults30_MIC[2] = 100 * abs(breakIndexes30MIC[1] - bAt30For_0_3n) / 30
  }
  
  if(breakIndexes30MIC[3] > 0){
    breakResults30_MIC[3] = 100 * abs(breakIndexes30MIC[3] - bAt30For_0_5n) / 30
  }
  
  if(breakIndexes30MIC[4] > 0){
    breakResults30_MIC[4] = 100 * abs(breakIndexes30MIC[4] - bAt30For_0_75n) / 30
  }
  
  # Calculating the error when n=100 (in terms of percent of total samples)
  if(breakIndexes100MIC[1] > 0){
    breakResults100_MIC[1] = 100 * abs(breakIndexes100MIC[1] - bAt100For_0_2n) / 100
  }
  
  if(breakIndexes100MIC[2] > 0){
    breakResults100_MIC[2] = 100 * abs(breakIndexes100MIC[2] - bAt100For_0_3n) / 100
  }
  
  if(breakIndexes100MIC[3] > 0){
    breakResults100_MIC[3] = 100 * abs(breakIndexes100MIC[3] - bAt100For_0_5n) / 100
  }
  
  if(breakIndexes100MIC[4] > 0){
    breakResults100_MIC[4] = 100 * abs(breakIndexes100MIC[4] - bAt100For_0_75n) / 100
  }
  
  # Calculating the error when n=200 (in terms of percent of total samples)
  if(breakIndexes200MIC[1] > 0){
    breakResults200_MIC[1] = 100 * abs(breakIndexes200MIC[1] - bAt200For_0_2n) / 200
  }
  
  if(breakIndexes200MIC[2] > 0){
    breakResults200_MIC[2] = 100 * abs(breakIndexes200MIC[2] - bAt200For_0_3n) / 200
  }
  
  if(breakIndexes200MIC[3] > 0){
    breakResults200_MIC[3] = 100 * abs(breakIndexes200MIC[3] - bAt200For_0_5n) / 200
  }
  
  if(breakIndexes200MIC[4] > 0){
    breakResults200_MIC[4] = 100 * abs(breakIndexes200MIC[4] - bAt200For_0_75n) / 200
  }
  
  breakResults.dataMIC <- data.frame(
    breakLocs = c("K=0.2n", "K=0.3n", "K=0.5n", "K=0.75n"),
    nEquals30 = breakResults30_MIC,
    nEquals100 = breakResults100_MIC,
    nEquals200 = breakResults200_MIC
  )
  
}
################################################################################
## Printing the output tables
if(runLR == 1){
  print('Results for the Likelihood ratio')
  print(results.dataLR)
  
  
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # Printing new line separator
    writeLines('\r\n')
    print("Percent error in detected break location for LR method")
    print(breakResults.dataLR)
  }
  
  # Printing new line separator
  writeLines('\r\n')
}
if(runSIC == 1){
  print("Coverage probabilities for the SIC method")
  print(results.dataSIC)
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # Printing new line separator
    writeLines('\r\n')
    print("Percent error in detected break location for SIC method")
    print(breakResults.dataSIC)
  }
}

if(runMIC == 1){
  print("Coverage probabilities for the MIC method")
  print(results.dataMIC)
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # Printing new line separator
    writeLines('\r\n')
    print("Percent error in detected break location for MIC method")
    print(breakResults.dataMIC)
  }
}


# Plotting power graphs for the techniques
if(showCoverageProbabilityCurves == 1){
  
  library(latex2exp)
  
  indexesToCalculate = seq(from=covProbParamStepSize, to=1-covProbParamStepSize, by=covProbParamStepSize)
  
  if(runLR==1){
    n30_PlotResultsLR <- matrix(0, 1,length(indexesToCalculate))
    n100_PlotResultsLR <- matrix(0, 1,length(indexesToCalculate))
    n200_PlotResultsLR <- matrix(0, 1,length(indexesToCalculate))
  }
  
  if(runSIC==1){
    n30_PlotResultsSIC <- matrix(0, 1, length(indexesToCalculate))
    n100_PlotResultsSIC <- matrix(0, 1, length(indexesToCalculate))
    n200_PlotResultsSIC <- matrix(0, 1, length(indexesToCalculate))
  }
  
  if(runMIC==1){
    n30_PlotResultsMIC <- matrix(0, 1, length(indexesToCalculate))
    n100_PlotResultsMIC <- matrix(0, 1, length(indexesToCalculate))
    n200_PlotResultsMIC <- matrix(0, 1, length(indexesToCalculate))
  }
  
  # Values for n = 100
  for (m in c(1:length(indexesToCalculate))){
    # Method Plot generation
    
    ## n=30
    simLen <- 30
    breakLocs <- indexesToCalculate[m] * simLen
    simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed, distFamily='E')
    if(runLR==1){
      # SIC Calc
      retVal = expIndepLRCalc(simMatrix, alpha=significanceLevel)
      # Grabbing the methods power
      n30_PlotResultsLR[1,m] <- retVal[1]
    }
    if(runSIC==1){
      retVal = expIndepSICCalc(simMatrix, alpha=significanceLevel)
      # Grabbing the coverage probability
      n30_PlotResultsSIC[1,m] = retVal[1]
    }
    if(runMIC==1){
      retVal = expIndepMICCalc(simMatrix, alpha=significanceLevel)
      # Grabbing the coverage probability
      n30_PlotResultsMIC[1,m] = retVal[1]
    }
    
    ## n=100
    simLen <- 100
    breakLocs <- indexesToCalculate[m] * simLen
    simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed, distFamily='E')
    if(runLR==1){
      # SIC Calc
      retVal = expIndepLRCalc(simMatrix, alpha=significanceLevel)
      # Grabbing the methods power
      n100_PlotResultsLR[1,m] <- retVal[1]
    }
    if(runSIC==1){
      retVal = expIndepSICCalc(simMatrix, alpha=significanceLevel)
      # Grabbing the coverage probability
      n100_PlotResultsSIC[1,m] = retVal[1]
    }
    if(runMIC==1){
      retVal = expIndepMICCalc(simMatrix, alpha=significanceLevel)
      # Grabbing the coverage probability
      n100_PlotResultsMIC[1,m] = retVal[1]
    }
    
    # Method Plot generation
    ## n=200
    simLen <- 200
    breakLocs <- indexesToCalculate[m] * simLen
    simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed, distFamily='E')
    if(runLR==1){
      # SIC Calc
      retVal = expIndepLRCalc(simMatrix, alpha=significanceLevel)
      # Grabbing the methods power
      n200_PlotResultsLR[1,m] <- retVal[1]
    }
    if(runSIC==1){
      retVal = expIndepSICCalc(simMatrix, alpha=significanceLevel)
      # Grabbing the coverage probability
      n200_PlotResultsSIC[1,m] = retVal[1]
    }
    if(runMIC==1){
      retVal = expIndepMICCalc(simMatrix, alpha=significanceLevel)
      # Grabbing the coverage probability
      n200_PlotResultsMIC[1,m] = retVal[1]
    }
  }
  # If all three are here, then plot on the same graph, otherwise plot on separate graphs
  if(runLR==1 && runSIC==1 && runMIC){
    ## Actual plotting code for n=30
    ylimits = c(0, max(c(n30_PlotResultsSIC[1,], n30_PlotResultsLR[1,])) + mean(c(n30_PlotResultsSIC[1,], n30_PlotResultsLR[1,])))
    plot(100*indexesToCalculate,n30_PlotResultsLR[1,],xlim=c(0,100), ylim=ylimits, col='blue', pch=20, main=TeX(paste('Power vs Break Index (n=30, $\\theta_1 = ', toString(param1), '$, $\\theta_2=', toString(param2), '$)')), ylab="Empirical Power", xlab="Break Index %")
    lines(100*indexesToCalculate,n30_PlotResultsLR[1,], col='blue')
    par(new=TRUE)
    plot(100*indexesToCalculate, n30_PlotResultsSIC[1,],xlim=c(0,100), ylim=ylimits,col='red',axes=FALSE, xlab='', ylab='', pch=20)
    lines(100*indexesToCalculate,n30_PlotResultsSIC[1,], col='red')
    par(new=TRUE)
    plot(100*indexesToCalculate, n30_PlotResultsMIC[1,],xlim=c(0,100), ylim=ylimits,col='green',axes=FALSE, xlab='', ylab='', pch=20)
    lines(100*indexesToCalculate,n30_PlotResultsMIC[1,], col='green')
    legend('topleft', c('Likelihood Ratio', 'SIC', 'MIC'), col=c('blue', 'red', 'green'), lty=c(1,1,1), cex=0.8)
    
    ## Actual plotting code for n=100
    ylimits = c(0, max(c(n100_PlotResultsSIC[1,], n100_PlotResultsLR[1,])) + mean(c(n100_PlotResultsSIC[1,], n100_PlotResultsLR[1,])))
    plot(100*indexesToCalculate,n100_PlotResultsLR[1,],xlim=c(0,100), ylim=c(0,1), col='blue', pch=20, main=TeX(paste('Power vs Break Index (n=100, $\\theta_1 = ', toString(param1), '$, $\\theta_2=', toString(param2), '$)')), ylab="Empirical Power", xlab="Break Index %")
    lines(100*indexesToCalculate,n100_PlotResultsLR[1,], col='blue')
    par(new=TRUE)
    plot(100*indexesToCalculate, n100_PlotResultsSIC[1,],xlim=c(0,100), ylim=c(0,1),col='red',axes=FALSE, xlab='', ylab='', pch=20)
    lines(100*indexesToCalculate, n100_PlotResultsSIC[1,], col='red')
    par(new=TRUE)
    plot(100*indexesToCalculate, n100_PlotResultsMIC[1,],xlim=c(0,100), ylim=c(0,1),col='green',axes=FALSE, xlab='', ylab='', pch=20)
    lines(100*indexesToCalculate, n100_PlotResultsMIC[1,], col='green')
    legend('topleft', c('Likelihood Ratio', 'SIC', 'MIC'), col=c('blue', 'red', 'green'), lty=c(1,1,1), cex=0.8)
    
    ## Actual plotting code for n=200
    ylimits = c(0, max(c(n200_PlotResultsSIC[1,], n200_PlotResultsLR[1,])) + mean(c(n200_PlotResultsSIC[1,], n200_PlotResultsLR[1,])))
    plot(100*indexesToCalculate,n200_PlotResultsLR[1,],xlim=c(0,100), ylim=c(0,1), col='blue', pch=20, main=TeX(paste('Power vs Break Index (n=200, $\\theta_1 = ', toString(param1), '$, $\\theta_2=', toString(param2), '$)')), ylab="Empirical Power", xlab="Break Index %")
    lines(100*indexesToCalculate,n200_PlotResultsLR[1,], col='blue')
    par(new=TRUE)
    plot(100*indexesToCalculate, n200_PlotResultsSIC[1,],xlim=c(0,100), ylim=c(0,1), col='red',axes=FALSE, xlab='', ylab='', pch=20)
    lines(100*indexesToCalculate,n200_PlotResultsSIC[1,], col='red')
    par(new=TRUE)
    plot(100*indexesToCalculate, n200_PlotResultsMIC[1,],xlim=c(0,100), ylim=c(0,1), col='green',axes=FALSE, xlab='', ylab='', pch=20)
    lines(100*indexesToCalculate,n200_PlotResultsMIC[1,], col='green')
    legend('topleft', c('Likelihood Ratio', 'SIC', 'MIC'), col=c('blue', 'red', 'green'), lty=c(1,1,1), cex=0.8)
  }
  else if(runLR==1){
    ## Plotting the the power versus break index, and connecting the points with lines (LR Method)
    plot(n100_PlotResultsLR[1,], main=TeX(paste('Power vs Break Index of LR Method (n=100, $\\theta_1 = ', toString(param1), '$, $\\theta_2=', toString(param2), '$)')), ylab="Empirical Power", xlab="Break Index %")
    lines(n100_PlotResultsLR[1,])
  }
  else if(runSIC==1){
    plot(n100_PlotResultsSIC[1,], main=TeX(paste('Power vs Break Index of SIC Method (n=100, $\\theta_1 = ', toString(param1), '$, $\\theta_2=', toString(param2), '$)')), ylab="Empirical Power", xlab="Break Index %")
    lines(n100_PlotResultsSIC[1,])
  }
  else if(runMIC==1){
    plot(n100_PlotResultsMIC[1,], main=TeX(paste('Power vs Break Index of MIC Method (n=100, $\\theta_1 = ', toString(param1), '$, $\\theta_2=', toString(param2), '$)')), ylab="Empirical Power", xlab="Break Index %")
    lines(n100_PlotResultsMIC[1,])
  }
}


# Code to implement the parameter sweep plot
if(showParamSweepPlot == 1){
  # Holding param1 constant, while we sweep out parameter 2
  pSwParam1 = 1
  
  # Only doing these simulations for n=100, and with k = 10, 20, and 50
  n = 100
  k = c(0.1*n, 0.2*n, 0.5*n)
  
  # Starting from 1/5 and going to 5
  # param2Vals = seq(from=0.2, to=5, by=paramSweepPlotStepSize)
  param2Vals = c(seq(from=0.2, to=1, length.out=round(paramSweepPlotStepNumber/2)), seq(from=1, to=5, length.out=round(paramSweepPlotStepNumber/2)))
  
  # Preallocating the appropriate simulation matrices
  if(runLR == 1){
    paramSweepResults.LR <- matrix(data=NA, nrow=length(param2Vals), ncol=length(k))
  }
  if(runSIC == 1){
    paramSweepResults.SIC <- matrix(data=NA, nrow=length(param2Vals), ncol=length(k))
  }
  if(runMIC == 1){
    paramSweepResults.MIC <- matrix(data=NA, nrow=length(param2Vals), ncol=length(k))
  }
  
  # Running a sweep for each of the three break locations
  for (j in 1:length(k)){
    # Sweeping over the param values
    for (m in 1:length(param2Vals)){
      simMatrix <-  indepDatagen(simLen=n, numSims=numSims, breakLocations=k[j], param1=pSwParam1, param2=param2Vals[m], seed=curSeed, distFamily='E')
      if(runLR == 1){
        paramSweepResults.LR[m,j] <- expIndepLRCalc(simMatrix, alpha=significanceLevel)[1]
      }
      if(runSIC == 1){
        paramSweepResults.SIC[m,j] <- expIndepSICCalc(simMatrix, alpha=significanceLevel)[1]
      }
      if(runMIC == 1){
        paramSweepResults.MIC[m,j] <- expIndepMICCalc(simMatrix, alpha=significanceLevel)[1]
      }
    }
  }

  # Start of plotting code
  if(runLR==1){
    
    plot(x=param2Vals, y=paramSweepResults.LR[,1], log='x', col=rgb(1,0,0), type='l', main='Likelihood Ratio method', xlab=TeX('$\\rho$'), ylab='Estimated Power')
    par(new=TRUE)
    plot(x=param2Vals, y=paramSweepResults.LR[,2], log='x', col=rgb(0,1,0), main='', xlab='', ylab='', axes=FALSE, type='l')
    par(new=TRUE)
    plot(x=param2Vals, y=paramSweepResults.LR[,3], log='x', col=rgb(0,0,1), main='', xlab='', ylab='', axes=FALSE, type='l')
    legend('top', c('k=0.1', 'k=0.2', 'k=0.5'), col=c(rgb(1,0,0), rgb(0,1,0), rgb(0,0,1)), lty=c(1,1,1), cex=0.8)
  }
  if(runSIC==1){
    plot(x=param2Vals, y=paramSweepResults.SIC[,1], log='x', col=rgb(1,0,0), type='l', main='SIC method', xlab=TeX('$\\rho$'), ylab='Estimated Power')
    par(new=TRUE)
    plot(x=param2Vals, y=paramSweepResults.SIC[,2], log='x', col=rgb(0,1,0), main='', xlab='', ylab='', axes=FALSE, type='l')
    par(new=TRUE)
    plot(x=param2Vals, y=paramSweepResults.SIC[,3], log='x', col=rgb(0,0,1), main='', xlab='', ylab='', axes=FALSE, type='l')
    legend('top', c('k=0.1', 'k=0.2', 'k=0.5'), col=c(rgb(1,0,0), rgb(0,1,0), rgb(0,0,1)), lty=c(1,1,1), cex=0.8)
  }
  if(runMIC==1){
    plot(x=param2Vals, y=paramSweepResults.MIC[,1], log='x', col=rgb(1,0,0), type='l', main='MIC method', xlab=TeX('$\\rho$'), ylab='Estimated Power')
    par(new=TRUE)
    plot(x=param2Vals, y=paramSweepResults.MIC[,2], log='x', col=rgb(0,1,0), main='', xlab='', ylab='', axes=FALSE, type='l')
    par(new=TRUE)
    plot(x=param2Vals, y=paramSweepResults.MIC[,3], log='x', col=rgb(0,0,1), main='', xlab='', ylab='', axes=FALSE, type='l')
    legend('top', c('k=0.1', 'k=0.2', 'k=0.5'), col=c(rgb(1,0,0), rgb(0,1,0), rgb(0,0,1)), lty=c(1,1,1), cex=0.8)
  }
}

