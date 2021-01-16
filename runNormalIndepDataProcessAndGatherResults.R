# Run normal independent data simulation and fill out results table

# Importing datagen function
source('indepDatagen.R')

# Importing CUSUM calculation function
source('CUSUM/CUSUMBreakEstimator.R')

# Initializing variables

## Arrays to store results of likelihood ratio calculation
nEquals50LR = c(0,0,0,0)
nEquals100LR = c(0,0,0,0)
nEquals200LR = c(0,0,0,0)

nEquals50CS = c(0,0,0,0)
nEquals100CS = c(0,0,0,0)
nEquals200CS = c(0,0,0,0)

## Parameters for the likelihood ratio simulation 
param1 = 0.1
param2 = 0.1
numSims = 100

# Performing simulations
################################################################################
## n = 50
simLen <- 50
### K=0.2n
breakLocs <- 0.2 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
# LR Calculation
source('LikelihoodRatio/normalDataParametricLikelihood.R')
nEquals50LR[1] = coverageProbability

# Cusum Calc
nEquals50CS[1] = CUSUMCalc(simMatrix)


### K=0.3n
breakLocs <- 0.3 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
# LR Calculation
source('LikelihoodRatio/normalDataParametricLikelihood.R')
nEquals50LR[2] = coverageProbability

# Cusum Calc
nEquals50CS[2] = CUSUMCalc(simMatrix)


### K=0.5n
breakLocs <- 0.5 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
# LR Calculation
source('LikelihoodRatio/normalDataParametricLikelihood.R')
nEquals50LR[3] = coverageProbability

# Cusum Calc
nEquals50CS[3] = CUSUMCalc(simMatrix)


### K=0.75n
breakLocs <- 0.75 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
# LR Calculation
source('LikelihoodRatio/normalDataParametricLikelihood.R')
nEquals50LR[4] = coverageProbability

# Cusum Calc
nEquals50CS[4] = CUSUMCalc(simMatrix)



################################################################################
## n = 100
simLen <- 100
### K=0.2n
breakLocs <- 0.2 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
# LR Calculation
source('LikelihoodRatio/normalDataParametricLikelihood.R')
nEquals100LR[1] = coverageProbability

# Cusum Calc
nEquals100CS[1] = CUSUMCalc(simMatrix)

### K=0.3n
breakLocs <- 0.3 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
# LR Calculation
source('LikelihoodRatio/normalDataParametricLikelihood.R')
nEquals100LR[2] = coverageProbability

# Cusum Calc
nEquals100CS[2] = CUSUMCalc(simMatrix)

### K=0.5n
breakLocs <- 0.5 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
# LR Calculation
source('LikelihoodRatio/normalDataParametricLikelihood.R')
nEquals100LR[3] = coverageProbability

# Cusum Calc
nEquals100CS[3] = CUSUMCalc(simMatrix)

### K=0.75n
breakLocs <- 0.75 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
# LR Calculation
source('LikelihoodRatio/normalDataParametricLikelihood.R')
nEquals100LR[4] = coverageProbability

# Cusum Calc
nEquals100CS[4] = CUSUMCalc(simMatrix)

################################################################################
## n = 200
simLen <- 200
### K=0.2n
breakLocs <- 0.2 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
# LR Calculation
source('LikelihoodRatio/normalDataParametricLikelihood.R')
nEquals200LR[1] = coverageProbability

# Cusum Calc
nEquals200CS[1] = CUSUMCalc(simMatrix)

### K=0.3n
breakLocs <- 0.3 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
# LR Calculation
source('LikelihoodRatio/normalDataParametricLikelihood.R')
nEquals200LR[2] = coverageProbability

# Cusum Calc
nEquals200CS[2] = CUSUMCalc(simMatrix)

### K=0.5n
breakLocs <- 0.5 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
# LR Calculation
source('LikelihoodRatio/normalDataParametricLikelihood.R')
nEquals200LR[3] = coverageProbability

# Cusum Calc
nEquals200CS[3] = CUSUMCalc(simMatrix)

### K=0.75n
breakLocs <- 0.75 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
# LR Calculation
source('LikelihoodRatio/normalDataParametricLikelihood.R')
nEquals200LR[4] = coverageProbability

# Cusum Calc
nEquals200CS[4] = CUSUMCalc(simMatrix)


################################################################################
# Constructing the output table for Likelihood Ratio
results.dataLR <- data.frame(
  breakLocs = c("K=0.2n", "K=0.3n", "K=0.5n", "K=0.75n"),
  nEquals50 = nEquals50LR,
  nEquals100 = nEquals100LR,
  nEquals200 = nEquals200LR
)

results.dataCS <- data.frame(
  breakLocs = c("K=0.2n", "K=0.3n", "K=0.5n", "K=0.75n"),
  nEquals50 = nEquals50CS,
  nEquals100 = nEquals100CS,
  nEquals200 = nEquals200CS
)


print('Results for the Likelihood ratio')
print(results.dataLR)

# Printing new line separator
writeLines('\r\n')

print("Results for the CUSUM calculation")
print(results.dataCS)


