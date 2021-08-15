
## This script generates multiple datasets of independent exponentially
## distributed data, with breaks at specified locations. The script
## computes the ability of two methods to detect the presence of these
## breaks: Conditional Likelihood Ratio, and Independent Likelihood Ratio.

## This code is effectively very similar to the runNormalDataProcessAndGatherResults
## code, but with the exception that we use exponentially distributed data 
## instead of normally distributed data

# Importing datagen function
source('timeSeriesDatagen.R')

# Importing the Likelihood Ration function
source('LikelihoodRatio/normalDataConditionalLikelihood.R')

# Importing ILR calculation function
source('LikelihoodRatio/normalDataParametricLikelihood.R')

# Initializing variables

# Modifiable Parameters
################################################################################
# Difference in mean between distributions
param1 = 0.5
param2 = 0.5
# Number of simulations to perform
numSims = 1000

# Seed for simulations (if set to 0, the seed is random, and any previously 
# set seed is globally cleared)
curSeed = 0

## Parameters for the Likelihood ratio, SIC and ILR algs.
significanceLevel = 0.05
# Whether to run LR calculation in this expirement or not
runLR = 1

# Whether to run ILR calculations in this experiment or not
runIndepLR = 1

## Parameters defining the type of output the user wants to see
# This will tell the code to either output the percent error in detected break locations, or not (1 to run, 0 to not)
showBreakLocDetectAcc = 0

################################################################################
## Arrays to store results of likelihood ratio and SIC calculations
nEquals30LR = c(0,0,0,0)
nEquals100LR = c(0,0,0,0)
nEquals200LR = c(0,0,0,0)

nEquals30ILR = c(0,0,0,0)
nEquals100ILR = c(0,0,0,0)
nEquals200ILR = c(0,0,0,0)

## Variables to store the average break location detected for each technique
if(showBreakLocDetectAcc){
  if(runLR){
    breakIndexes30LR = c(0, 0, 0, 0)
    breakIndexes100LR = c(0, 0, 0, 0)
    breakIndexes200LR = c(0, 0, 0, 0) 
  }
  if(runIndepLR){
    breakIndexes30ILR = c(0,0,0,0)
    breakIndexes100ILR = c(0,0,0,0)
    breakIndexes200ILR = c(0,0,0,0)
  }
}
# Performing simulations
################################################################################
## n = 30
simLen <- 30
### K=0.2n
breakLocs <- 0.2 * simLen
simMatrix <- timeSeriesDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation
if(runLR==1){
  retVal = normalCondLRCalc(simMatrix, alpha=significanceLevel, mu=0, sigma=1)
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

if(runIndepLR==1){
  # ILR Calc
  retVal = normalIndepLRCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals30ILR[1] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes30ILR[1] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes30ILR[1] = -1
    }
  }
}


### K=0.3n
breakLocs <- 0.3 * simLen
simMatrix <- timeSeriesDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation

if(runLR==1){
  retVal = normalCondLRCalc(simMatrix, alpha=significanceLevel, mu=0, sigma=1)
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

if(runIndepLR==1){
  # SIC Calc
  retVal = normalIndepLRCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals30ILR[2] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes30ILR[2] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes30ILR[2] = -1
    }
  }
}

### K=0.5n
breakLocs <- 0.5 * simLen
simMatrix <- timeSeriesDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)

# LR Calculation
if(runLR==1){
  retVal = normalCondLRCalc(simMatrix, alpha=significanceLevel, mu=0, sigma=1)
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

if(runIndepLR==1){
  # ILR Calc
  retVal = normalIndepLRCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals30ILR[3] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes30ILR[3] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes30ILR[3] = -1
    }
  }
}

### K=0.75n
breakLocs <- 0.75 * simLen
simMatrix <- timeSeriesDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation
if(runLR==1){
  retVal = normalCondLRCalc(simMatrix, alpha=significanceLevel, mu=0, sigma=1)
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

if(runIndepLR==1){
  # ILR Calc
  retVal = normalIndepLRCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals30ILR[4] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes30ILR[4] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes30ILR[4] = -1
    }
  }
}



################################################################################
## n = 100
simLen <- 100
### K=0.2n
breakLocs <- 0.2 * simLen
simMatrix <- timeSeriesDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation
if(runLR==1){
  retVal = normalCondLRCalc(simMatrix, alpha=significanceLevel, mu=0, sigma=1)
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

if(runIndepLR==1){
  # ILR Calc
  retVal = normalIndepLRCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals100ILR[1] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes100ILR[1] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes100ILR[1] = -1
    }
  }
}


### K=0.3n
breakLocs <- 0.3 * simLen
simMatrix <- timeSeriesDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation

if(runLR==1){
  retVal = normalCondLRCalc(simMatrix, alpha=significanceLevel, mu=0, sigma=1)
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

if(runIndepLR==1){
  # SIC Calc
  retVal = normalIndepLRCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals100ILR[2] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes100ILR[2] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes100ILR[2] = -1
    }
  }
}

### K=0.5n
breakLocs <- 0.5 * simLen
simMatrix <- timeSeriesDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)

# LR Calculation
if(runLR==1){
  retVal = normalCondLRCalc(simMatrix, alpha=significanceLevel, mu=0, sigma=1)
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

if(runIndepLR==1){
  # ILR Calc
  retVal = normalIndepLRCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals100ILR[3] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes100ILR[3] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes100ILR[3] = -1
    }
  }
}

### K=0.75n
breakLocs <- 0.75 * simLen
simMatrix <- timeSeriesDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation
if(runLR==1){
  retVal = normalCondLRCalc(simMatrix, alpha=significanceLevel, mu=0, sigma=1)
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

if(runIndepLR==1){
  # ILR Calc
  retVal = normalIndepLRCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals100ILR[4] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes100ILR[4] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes100ILR[4] = -1
    }
  }
}

## n = 200
simLen <- 200
### K=0.2n
breakLocs <- 0.2 * simLen
simMatrix <- timeSeriesDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation
if(runLR==1){
  retVal = normalCondLRCalc(simMatrix, alpha=significanceLevel, mu=0, sigma=1)
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

if(runIndepLR==1){
  # ILR Calc
  retVal = normalIndepLRCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals200ILR[1] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes200ILR[1] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes200ILR[1] = -1
    }
  }
}


### K=0.3n
breakLocs <- 0.3 * simLen
simMatrix <- timeSeriesDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation

if(runLR==1){
  retVal = normalCondLRCalc(simMatrix, alpha=significanceLevel, mu=0, sigma=1)
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

if(runIndepLR==1){
  # SIC Calc
  retVal = normalIndepLRCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals200ILR[2] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes200ILR[2] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes200ILR[2] = -1
    }
  }
}

### K=0.5n
breakLocs <- 0.5 * simLen
simMatrix <- timeSeriesDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)

# LR Calculation
if(runLR==1){
  retVal = normalCondLRCalc(simMatrix, alpha=significanceLevel, mu=0, sigma=1)
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

if(runIndepLR==1){
  # ILR Calc
  retVal = normalIndepLRCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals200ILR[3] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes200ILR[3] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes200ILR[3] = -1
    }
  }
}

### K=0.75n
breakLocs <- 0.75 * simLen
simMatrix <- timeSeriesDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2, seed=curSeed)
# LR Calculation
if(runLR==1){
  retVal = normalCondLRCalc(simMatrix, alpha=significanceLevel, mu=0, sigma=1)
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

if(runIndepLR==1){
  # ILR Calc
  retVal = normalIndepLRCalc(simMatrix, alpha=significanceLevel)
  # Grabbing the coverage probability
  nEquals200ILR[4] = retVal[1]
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # If we detected at least one break then we can get an estimate of the break index
    if(retVal[1] > 0){
      # Calculating the average index of detected break
      indexesOfDetectedBreaks = retVal[-1]
      # Removing indexes where a break wasn't detected
      indexesOfDetectedBreaks = indexesOfDetectedBreaks[indexesOfDetectedBreaks != -1]
      breakIndexes200ILR[4] = mean(indexesOfDetectedBreaks)
    }else{
      breakIndexes200ILR[4] = -1
    }
  }
}


################################################################################
## Constructing the output table for Likelihood Ratio and SIC
results.dataCLR <- data.frame(
  breakLocs = c("K=0.2n", "K=0.3n", "K=0.5n", "K=0.75n"),
  nEquals30 = 1 - nEquals30LR,
  nEquals100 = 1 - nEquals100LR,
  nEquals200 = 1 - nEquals200LR
)


results.dataILR <- data.frame(
  breakLocs = c("K=0.2n", "K=0.3n", "K=0.5n", "K=0.75n"),
  nEquals30 = 1 - nEquals30ILR,
  nEquals100 = 1 - nEquals100ILR,
  nEquals200 = 1 - nEquals200ILR
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
  
  breakResults.dataCLR <- data.frame(
    breakLocs = c("K=0.2n", "K=0.3n", "K=0.5n", "K=0.75n"),
    nEquals30 = 1 - breakResults30_LR,
    nEquals100 = 1 - breakResults100_LR,
    nEquals200 = 1 - breakResults200_LR
  )
  
  ## Calculating the error of the ILR's break location detection
  breakResults30_ILR = c(-1,-1,-1,-1)
  breakResults100_ILR = c(-1, -1, -1, -1)
  breakResults200_ILR = c(-1,-1,-1,-1)
  # Calculating the error when n=30 (in terms of percent of total samples)
  # If a break was detected we can do the calculation, otherwise we leave it as -1
  if(breakIndexes30ILR[1] > 0){
    breakResults30_ILR[1] = 100 * abs(breakIndexes30ILR[1] - bAt30For_0_2n) / 30
  }
  
  if(breakIndexes30ILR[2] > 0){
    breakResults30_ILR[2] = 100 * abs(breakIndexes30ILR[1] - bAt30For_0_3n) / 30
  }
  
  if(breakIndexes30ILR[3] > 0){
    breakResults30_ILR[3] = 100 * abs(breakIndexes30ILR[3] - bAt30For_0_5n) / 30
  }
  
  if(breakIndexes30ILR[4] > 0){
    breakResults30_ILR[4] = 100 * abs(breakIndexes30ILR[4] - bAt30For_0_75n) / 30
  }
  
  # Calculating the error when n=100 (in terms of percent of total samples)
  if(breakIndexes100ILR[1] > 0){
    breakResults100_ILR[1] = 100 * abs(breakIndexes100ILR[1] - bAt100For_0_2n) / 100
  }
  
  if(breakIndexes100ILR[2] > 0){
    breakResults100_ILR[2] = 100 * abs(breakIndexes100ILR[2] - bAt100For_0_3n) / 100
  }
  
  if(breakIndexes100ILR[3] > 0){
    breakResults100_ILR[3] = 100 * abs(breakIndexes100ILR[3] - bAt100For_0_5n) / 100
  }
  
  if(breakIndexes100ILR[4] > 0){
    breakResults100_ILR[4] = 100 * abs(breakIndexes100ILR[4] - bAt100For_0_75n) / 100
  }
  
  # Calculating the error when n=200 (in terms of percent of total samples)
  if(breakIndexes200ILR[1] > 0){
    breakResults200_ILR[1] = 100 * abs(breakIndexes200ILR[1] - bAt200For_0_2n) / 200
  }
  
  if(breakIndexes200ILR[2] > 0){
    breakResults200_ILR[2] = 100 * abs(breakIndexes200ILR[2] - bAt200For_0_3n) / 200
  }
  
  if(breakIndexes200ILR[3] > 0){
    breakResults200_ILR[3] = 100 * abs(breakIndexes200ILR[3] - bAt200For_0_5n) / 200
  }
  
  if(breakIndexes200ILR[4] > 0){
    breakResults200_ILR[4] = 100 * abs(breakIndexes200ILR[4] - bAt200For_0_75n) / 200
  }
  
  breakResults.dataILR <- data.frame(
    breakLocs = c("K=0.2n", "K=0.3n", "K=0.5n", "K=0.75n"),
    nEquals30 = breakResults30_ILR,
    nEquals100 = breakResults100_ILR,
    nEquals200 = breakResults200_ILR
  )
  
  
}

################################################################################
## Printing the output tables
if(runLR == 1){
  print('Results for the Conditional LR method')
  print(results.dataCLR)
  
  
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # Printing new line separator
    writeLines('\r\n')
    print("Percent error in detected break location for CLR method")
    print(breakResults.dataCLR)
  }
  
  # Printing new line separator
  writeLines('\r\n')
}

if(runIndepLR == 1){
  print("Coverage probabilities for the Independent LR method")
  print(results.dataILR)
  # Only running if we want to calculate the break location detection accuracy
  if(showBreakLocDetectAcc == 1){
    # Printing new line separator
    writeLines('\r\n')
    print("Percent error in detected break location for ILR method")
    print(breakResults.dataILR)
  }
}


