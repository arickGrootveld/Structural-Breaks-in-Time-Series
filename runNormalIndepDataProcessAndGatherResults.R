
## This script generates multiple datasets of independent normally 
## distributed data, with breaks at specified locations. The script
## computes the ability of two methods to detect the pressence of these
## breaks: Likelihood Ratio, and CUSUM.

# Importing datagen function
source('indepDatagen.R')

# Importing the Likelihood Ration function
source('LikelihoodRatio/normalDataParametricLikelihood.R')

# Importing CUSUM calculation function
source('CUSUM/CUSUMBreakEstimator.R')

# Initializing variables

# Modifiable Parameters
################################################################################
# Difference in mean between distributions
param1 = 1.5
param2 = 2.5
# Number of simulations to perform
numSims = 5000

# Seed for simulations (if set to 0, the seed is random, and any previously 
# set seed is globally cleared)
curSeed = 0

## Parameters for the Likelihood ratio
# lr params set to the simulation params to improve accuracy, but can be adjusted
# to evaluate the amount of error caused by parametric mismatch
lrParam1 = param1
lrParam2 = param2
significanceLevel = 0.05
# Whether to run LR calculation in this expirement or not
runLR = 1

## Parameters for the CUSUM
criticalValueCS = 2.408
longRunVariance = 1
# Whether to run CUSUM calculations in this expirement or not
runCUSUM = 1

## Parameters defining the type of output the user wants to see
# This will tell the code to either output the percent error in detected break locations, or not (1 to run, 0 to not)
showBreakLocDetectAcc = 1

# This will tell the code to output coverage probability curves for the method (TODO: Implement this)
showCoverageProbabilityCurves = 0
# If the above variable is 1, then this variable will set the step size of the variable used to linearly interpolate the 
# line of the coverage probability curve(as a fraction percent out of 1, i.e.
# 0.1 would plot the score in intervals of 10% of the total sequence length)
covProbParamStepSize = 0.04

# Variable declarations
################################################################################
## Arrays to store results of likelihood ratio and CUSUM calculations
nEquals25LR = c(0,0,0,0)
nEquals100LR = c(0,0,0,0)
nEquals200LR = c(0,0,0,0)

nEquals25CS = c(0,0,0,0)
nEquals100CS = c(0,0,0,0)
nEquals200CS = c(0,0,0,0)

## Variables to store the average break location detected for each technique
breakIndexes25LR = c(0, 0, 0, 0)
breakIndexes100LR = c(0, 0, 0, 0)
breakIndexes200LR = c(0, 0, 0, 0) 

breakIndexes25CS = c(0,0,0,0)
breakIndexes100CS = c(0,0,0,0)
breakIndexes200CS = c(0,0,0,0)


# Performing simulations
################################################################################
## n = 25
simLen <- 25
### K=0.2n
breakLocs <- 0.2 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation
if(runLR==1){
  retVal = normalIndepLRCalc(lrParam1, lrParam2, simMatrix, alpha=significanceLevel)
  nEquals25LR[1] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes25LR[1] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes25LR[1] = -1
    }
  }
}

if(runCUSUM==1){
  # Cusum Calc
  retVal = CUSUMCalc(simMatrix, critVal=criticalValueCS, longRunVar=longRunVariance)
  # Grabbing the coverage probability
  nEquals25CS[1] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes25CS[1] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes25CS[1] = -1
    }
  }
}


### K=0.3n
breakLocs <- 0.3 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation

if(runLR==1){
  retVal = normalIndepLRCalc(lrParam1, lrParam2, simMatrix, alpha=significanceLevel)
  nEquals25LR[2] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes25LR[2] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes25LR[2] = -1
    }
  }
}

if(runCUSUM==1){
  # Cusum Calc
  retVal = CUSUMCalc(simMatrix, critVal=criticalValueCS, longRunVar=longRunVariance)
  # Grabbing the coverage probability
  nEquals25CS[2] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes25CS[2] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes25CS[2] = -1
    }
  }
}


### K=0.5n
breakLocs <- 0.5 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation

if(runLR==1){
  retVal = normalIndepLRCalc(lrParam1, lrParam2, simMatrix, alpha=significanceLevel)
  nEquals25LR[3] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes25LR[3] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes25LR[3] = -1
    }
  }
}

if(runCUSUM==1){
  # Cusum Calc
  retVal = CUSUMCalc(simMatrix, critVal=criticalValueCS, longRunVar=longRunVariance)
  # Grabbing the coverage probability
  nEquals25CS[3] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes25CS[3] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes25CS[3] = -1
    }
  }
}


### K=0.75n
breakLocs <- 0.75 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation
if(runLR==1){
  retVal = normalIndepLRCalc(lrParam1, lrParam2, simMatrix, alpha=significanceLevel)
  nEquals25LR[4] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes25LR[4] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes25LR[4] = -1
    }
  }
}

if(runCUSUM==1){
  # Cusum Calc
  retVal = CUSUMCalc(simMatrix, critVal=criticalValueCS, longRunVar=longRunVariance)
  # Grabbing the coverage probability
  nEquals25CS[4] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Getting all break detection indexes so we can calculate the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes25CS[4] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes25CS[4] = -1
    }
  }
}


################################################################################
## n = 100
simLen <- 100
### K=0.2n
breakLocs <- 0.2 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation

if(runLR==1){
  retVal = normalIndepLRCalc(lrParam1, lrParam2, simMatrix, alpha=significanceLevel)
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

if(runCUSUM==1){
  # Cusum Calc
  retVal = CUSUMCalc(simMatrix, critVal=criticalValueCS, longRunVar=longRunVariance)
  # Grabbing the coverage probability
  nEquals100CS[1] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes100CS[1] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes100CS[1] = -1
    }
  }
}

### K=0.3n
breakLocs <- 0.3 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation

if(runLR==1){
  retVal = normalIndepLRCalc(lrParam1, lrParam2, simMatrix, alpha=significanceLevel)
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

if(runCUSUM==1){
  # Cusum Calc
  retVal = CUSUMCalc(simMatrix, critVal=criticalValueCS, longRunVar=longRunVariance)
  # Grabbing the coverage probability
  nEquals100CS[2] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes100CS[2] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes100CS[2] = -1
    }
  }
}

### K=0.5n
breakLocs <- 0.5 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation

if(runLR==1){
  retVal = normalIndepLRCalc(lrParam1, lrParam2, simMatrix, alpha=significanceLevel)
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

if(runCUSUM==1){
  # Cusum Calc
  retVal = CUSUMCalc(simMatrix, critVal=criticalValueCS, longRunVar=longRunVariance)
  # Grabbing the coverage probability
  nEquals100CS[3] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes100CS[3] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes100CS[3] = -1
    }
  }
}

### K=0.75n
breakLocs <- 0.75 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation

if(runLR==1){
  retVal = normalIndepLRCalc(lrParam1, lrParam2, simMatrix, alpha=significanceLevel)
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

if(runCUSUM==1){
  # Cusum Calc
  retVal = CUSUMCalc(simMatrix, critVal=criticalValueCS, longRunVar=longRunVariance)
  # Grabbing the coverage probability
  nEquals100CS[4] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes100CS[4] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes100CS[4] = -1
    }
  }
}

################################################################################
## n = 200
simLen <- 200
### K=0.2n
breakLocs <- 0.2 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation

if(runLR==1){
  retVal = normalIndepLRCalc(lrParam1, lrParam2, simMatrix, alpha=significanceLevel)
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

if(runCUSUM==1){
  # Cusum Calc
  retVal = CUSUMCalc(simMatrix, critVal=criticalValueCS, longRunVar=longRunVariance)
  # Grabbing the coverage probability
  nEquals200CS[1] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes200CS[1] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes200CS[1] = -1
    }
  }
}
### K=0.3n
breakLocs <- 0.3 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation

if(runLR==1){
  retVal = normalIndepLRCalc(lrParam1, lrParam2, simMatrix, alpha=significanceLevel)
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

if(runCUSUM==1){
  # Cusum Calc
  retVal = CUSUMCalc(simMatrix, critVal=criticalValueCS, longRunVar=longRunVariance)
  # Grabbing the coverage probability
  nEquals200CS[2] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes200CS[2] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes200CS[2] = -1
    }
  }
}

### K=0.5n
breakLocs <- 0.5 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation

if(runLR==1){
  retVal = normalIndepLRCalc(lrParam1, lrParam2, simMatrix, alpha=significanceLevel)
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


if(runCUSUM==1){
  # Cusum Calc
  retVal = CUSUMCalc(simMatrix, critVal=criticalValueCS, longRunVar=longRunVariance)
  # Grabbing the coverage probability
  nEquals200CS[3] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes200CS[3] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes200CS[3] = -1
    }
  }
}

### K=0.75n
breakLocs <- 0.75 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation

if(runLR==1){
  retVal = normalIndepLRCalc(lrParam1, lrParam2, simMatrix, alpha=significanceLevel)
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


if(runCUSUM==1){
  # Cusum Calc
  retVal = CUSUMCalc(simMatrix, critVal=criticalValueCS, longRunVar=longRunVariance)
  # Grabbing the coverage probability
  nEquals200CS[4] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes200CS[4] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes200CS[4] = -1
    }
  }
}


################################################################################
## Constructing the output table for Likelihood Ratio and CUSUM
results.dataLR <- data.frame(
  breakLocs = c("K=0.2n", "K=0.3n", "K=0.5n", "K=0.75n"),
  nEquals25 = nEquals25LR,
  nEquals100 = nEquals100LR,
  nEquals200 = nEquals200LR
)

results.dataCS <- data.frame(
  breakLocs = c("K=0.2n", "K=0.3n", "K=0.5n", "K=0.75n"),
  nEquals25 = nEquals25CS,
  nEquals100 = nEquals100CS,
  nEquals200 = nEquals200CS
)

## Constructing the output table for the detected break locations
# Getting actual break locations so we can calculate error of the detected break locations
if(showBreakLocDetectAcc == 1){
  bAt25For_0_2n = round(25 * 0.2)
  bAt25For_0_3n = round(25 * 0.3)
  bAt25For_0_5n = round(25 * 0.5)
  bAt25For_0_75n = round(25 * 0.75)
  
  bAt100For_0_2n = round(100 * 0.2)
  bAt100For_0_3n = round(100 * 0.3)
  bAt100For_0_5n = round(100 * 0.5)
  bAt100For_0_75n = round(100 * 0.75)
  
  bAt200For_0_2n = round(200 * 0.2)
  bAt200For_0_3n = round(200 * 0.3)
  bAt200For_0_5n = round(200 * 0.5)
  bAt200For_0_75n = round(200 * 0.75)
  
  
  ## Calculating the error of the LR method's break location detection
  breakResults25_LR = c(-1,-1,-1,-1)
  breakResults100_LR = c(-1, -1, -1, -1)
  breakResults200_LR = c(-1,-1,-1,-1)
  # Calculating the error when n=25 (in terms of percent of total samples)
  # If a break was detected we can do the calculation, otherwise we leave it as -1
  if(breakIndexes25LR[1] > 0){
    breakResults25_LR[1] = 100 * abs(breakIndexes25LR[1] - bAt25For_0_2n) / 25
  }
  
  if(breakIndexes25LR[2] > 0){
    breakResults25_LR[2] = 100 * abs(breakIndexes25LR[1] - bAt25For_0_3n) / 25
  }
  
  if(breakIndexes25LR[3] > 0){
    breakResults25_LR[3] = 100 * abs(breakIndexes25LR[3] - bAt25For_0_5n) / 25
  }
  
  if(breakIndexes25LR[4] > 0){
    breakResults25_LR[4] = 100 * abs(breakIndexes25LR[4] - bAt25For_0_75n) / 25
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
    nEquals25 = breakResults25_LR,
    nEquals100 = breakResults100_LR,
    nEquals200 = breakResults200_LR
  )
  
  
  ## Calculating the error of the CUSUM's break location detection
  breakResults25_CS = c(-1,-1,-1,-1)
  breakResults100_CS = c(-1, -1, -1, -1)
  breakResults200_CS = c(-1,-1,-1,-1)
  # Calculating the error when n=25 (in terms of percent of total samples)
  # If a break was detected we can do the calculation, otherwise we leave it as -1
  if(breakIndexes25CS[1] > 0){
    breakResults25_CS[1] = 100 * abs(breakIndexes25CS[1] - bAt25For_0_2n) / 25
  }
  
  if(breakIndexes25CS[2] > 0){
    breakResults25_CS[2] = 100 * abs(breakIndexes25CS[1] - bAt25For_0_3n) / 25
  }
  
  if(breakIndexes25CS[3] > 0){
    breakResults25_CS[3] = 100 * abs(breakIndexes25CS[3] - bAt25For_0_5n) / 25
  }
  
  if(breakIndexes25CS[4] > 0){
    breakResults25_CS[4] = 100 * abs(breakIndexes25CS[4] - bAt25For_0_75n) / 25
  }
  
  # Calculating the error when n=100 (in terms of percent of total samples)
  if(breakIndexes100CS[1] > 0){
    breakResults100_CS[1] = 100 * abs(breakIndexes100CS[1] - bAt100For_0_2n) / 100
  }
  
  if(breakIndexes100CS[2] > 0){
    breakResults100_CS[2] = 100 * abs(breakIndexes100CS[2] - bAt100For_0_3n) / 100
  }
  
  if(breakIndexes100CS[3] > 0){
    breakResults100_CS[3] = 100 * abs(breakIndexes100CS[3] - bAt100For_0_5n) / 100
  }
  
  if(breakIndexes100CS[4] > 0){
    breakResults100_CS[4] = 100 * abs(breakIndexes100CS[4] - bAt100For_0_75n) / 100
  }
  
  # Calculating the error when n=200 (in terms of percent of total samples)
  if(breakIndexes200CS[1] > 0){
    breakResults200_CS[1] = 100 * abs(breakIndexes200CS[1] - bAt200For_0_2n) / 200
  }
  
  if(breakIndexes200CS[2] > 0){
    breakResults200_CS[2] = 100 * abs(breakIndexes200CS[2] - bAt200For_0_3n) / 200
  }
  
  if(breakIndexes200CS[3] > 0){
    breakResults200_CS[3] = 100 * abs(breakIndexes200CS[3] - bAt200For_0_5n) / 200
  }
  
  if(breakIndexes200CS[4] > 0){
    breakResults200_CS[4] = 100 * abs(breakIndexes200CS[4] - bAt200For_0_75n) / 200
  }
  
  breakResults.dataCS <- data.frame(
    breakLocs = c("K=0.2n", "K=0.3n", "K=0.5n", "K=0.75n"),
    nEquals25 = breakResults25_CS,
    nEquals100 = breakResults100_CS,
    nEquals200 = breakResults200_CS
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
if(runCUSUM == 1){
  print("Coverage probabilities for the CUSUM method")
  print(results.dataCS)
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # Printing new line separator
    writeLines('\r\n')
    print("Percent error in detected break location for CUSUM method")
    print(breakResults.dataCS)
  }
}
# Plotting power graphs for the techniques
if(showCoverageProbabilityCurves == 1){
  
  library(latex2exp)
  
  indexesToCalculate = seq(from=covProbParamStepSize, to=1-covProbParamStepSize, by=covProbParamStepSize)
  
  if(runLR==1){
    n25_PlotResultsLR <- matrix(0, 1,length(indexesToCalculate))
    n100_PlotResultsLR <- matrix(0, 1,length(indexesToCalculate))
    n200_PlotResultsLR <- matrix(0, 1,length(indexesToCalculate))
  }
  
  if(runCUSUM==1){
    n25_PlotResultsCS <- matrix(0, 1, length(indexesToCalculate))
    n100_PlotResultsCS <- matrix(0, 1, length(indexesToCalculate))
    n200_PlotResultsCS <- matrix(0, 1, length(indexesToCalculate))
  }
  
  # Values for n = 100
  for (m in c(1:length(indexesToCalculate))){
    # Method Plot generation
    
    ## n=25
    simLen <- 25
    breakLocs <- indexesToCalculate[m] * simLen
    simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
    if(runLR==1){
      # Cusum Calc
      retVal = normalIndepLRCalc(lrParam1, lrParam2, simMatrix, alpha=significanceLevel)
      # Grabbing the methods power
      n25_PlotResultsLR[1,m] <- retVal[1]
    }
    if(runCUSUM==1){
      retVal = CUSUMCalc(simMatrix, critVal=criticalValueCS, longRunVar=longRunVariance)
      # Grabbing the coverage probability
      n25_PlotResultsCS[1,m] = retVal[1]
    }
    
    ## n=100
    simLen <- 100
    breakLocs <- indexesToCalculate[m] * simLen
    simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
    if(runLR==1){
      # Cusum Calc
      retVal = normalIndepLRCalc(lrParam1, lrParam2, simMatrix, alpha=significanceLevel)
      # Grabbing the methods power
      n100_PlotResultsLR[1,m] <- retVal[1]
    }
    if(runCUSUM==1){
      retVal = CUSUMCalc(simMatrix, critVal=criticalValueCS, longRunVar=longRunVariance)
      # Grabbing the coverage probability
      n100_PlotResultsCS[1,m] = retVal[1]
    }
    
    # Method Plot generation
    ## n=200
    simLen <- 200
    breakLocs <- indexesToCalculate[m] * simLen
    simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
    if(runLR==1){
      # Cusum Calc
      retVal = normalIndepLRCalc(lrParam1, lrParam2, simMatrix, alpha=significanceLevel)
      # Grabbing the methods power
      n200_PlotResultsLR[1,m] <- retVal[1]
    }
    if(runCUSUM==1){
      retVal = CUSUMCalc(simMatrix, critVal=criticalValueCS, longRunVar=longRunVariance)
      # Grabbing the coverage probability
      n200_PlotResultsCS[1,m] = retVal[1]
    }
  }
  # If both are here, then plot on the same graph, otherwise plot on seperate graphs
  if(runLR==1 && runCUSUM==1){
    ## Actual plotting code for n=25
    ylimits = c(0, max(c(n25_PlotResultsCS[1,], n25_PlotResultsLR[1,])) + mean(c(n25_PlotResultsCS[1,], n25_PlotResultsLR[1,])))
    plot(100*indexesToCalculate,n25_PlotResultsLR[1,],xlim=c(0,100), ylim=ylimits, col='blue', pch=20, main=TeX('Power vs Break Index (n=25, $\\theta_1 = 1.5$, $\\theta_2=2.5$)'), ylab="Empirical Power", xlab="Break Index %")
    lines(100*indexesToCalculate,n25_PlotResultsLR[1,], col='blue')
    par(new=TRUE)
    plot(100*indexesToCalculate, n25_PlotResultsCS[1,],xlim=c(0,100), ylim=ylimits,col='red',axes=FALSE, xlab='', ylab='', pch=20)
    lines(100*indexesToCalculate,n25_PlotResultsCS[1,], col='red')
    legend('topleft', c('Likelihood Ratio', 'CUSUM'), col=c('blue', 'red'), lty=c(1,1), cex=0.8)
    
    ## Actual plotting code for n=100
    ylimits = c(0, max(c(n100_PlotResultsCS[1,], n100_PlotResultsLR[1,])) + mean(c(n100_PlotResultsCS[1,], n100_PlotResultsLR[1,])))
    plot(100*indexesToCalculate,n100_PlotResultsLR[1,],xlim=c(0,100), ylim=c(0,1), col='blue', pch=20, main=TeX('Power vs Break Index (n=100, $\\theta_1 = 1.5$, $\\theta_2=2.5$)'), ylab="Empirical Power", xlab="Break Index %")
    lines(100*indexesToCalculate,n100_PlotResultsLR[1,], col='blue')
    par(new=TRUE)
    plot(100*indexesToCalculate, n100_PlotResultsCS[1,],xlim=c(0,100), ylim=c(0,1),col='red',axes=FALSE, xlab='', ylab='', pch=20)
    lines(100*indexesToCalculate, n100_PlotResultsCS[1,], col='red')
    legend('topleft', c('Likelihood Ratio', 'CUSUM'), col=c('blue', 'red'), lty=c(1,1), cex=0.8)
    
    ## Actual plotting code for n=200
    ylimits = c(0, max(c(n200_PlotResultsCS[1,], n200_PlotResultsLR[1,])) + mean(c(n200_PlotResultsCS[1,], n200_PlotResultsLR[1,])))
    plot(100*indexesToCalculate,n200_PlotResultsLR[1,],xlim=c(0,100), ylim=c(0,1), col='blue', pch=20, main=TeX('Power vs Break Index (n=200, $\\theta_1 = 1.5$, $\\theta_2=2.5$)'), ylab="Empirical Power", xlab="Break Index %")
    lines(100*indexesToCalculate,n200_PlotResultsLR[1,], col='blue')
    par(new=TRUE)
    plot(100*indexesToCalculate, n200_PlotResultsCS[1,],xlim=c(0,100), ylim=c(0,1), col='red',axes=FALSE, xlab='', ylab='', pch=20)
    lines(100*indexesToCalculate,n200_PlotResultsCS[1,], col='red')
    legend('topleft', c('Likelihood Ratio', 'CUSUM'), col=c('blue', 'red'), lty=c(1,1), cex=0.8)
  }
  else if(runLR==1){
    ## Plotting the the power versus break index, and connecting the points with lines (LR Method)
    plot(n100_PlotResultsLR[1,], main=TeX('Power vs Break Index of LR Method (n=100, $\\theta_1 = 1.5$, $\\theta_2=2.5$)'), ylab="Empirical Power", xlab="Break Index %")
    lines(n100_PlotResultsLR[1,])
  }
  else if(runCUSUM==1){
    plot(n100_PlotResultsCS[1,], main=TeX('Power vs Break Index of CUSUM Method (n=100, $\\theta_1 = 1.5$, $\\theta_2=2.5$)'), ylab="Empirical Power", xlab="Break Index %")
    lines(n100_PlotResultsCS[1,])
  }
  
  
  
  
}