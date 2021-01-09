# Run normal independent data simulation and fill out results table

# Importing datagen function
source('indepDatagen.R')

# Initializing variables

## Arrays to store results
nEquals50 = c(0,0,0,0)
nEquals100 = c(0,0,0,0)
nEquals200 = c(0,0,0,0)

## Parameters for the simulation
param1 = 0.1
param2 = 0.5
numSims = 10

# Performing simulations
################################################################################
## n = 50
simLen <- 50
### K=0.2n
breakLocs <- 0.2 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
source('normalDataParametricLikelihood.R')

nEquals50[1] = coverageProbability


### K=0.3n
breakLocs <- 0.3 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
source('normalDataParametricLikelihood.R')

nEquals50[2] = coverageProbability


### K=0.5n
breakLocs <- 0.5 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
source('normalDataParametricLikelihood.R')

nEquals50[3] = coverageProbability


### K=0.75n
breakLocs <- 0.75 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
source('normalDataParametricLikelihood.R')

nEquals50[4] = coverageProbability



################################################################################
## n = 100
simLen <- 100
### K=0.2n
breakLocs <- 0.2 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
source('normalDataParametricLikelihood.R')

nEquals100[1] = coverageProbability


### K=0.3n
breakLocs <- 0.3 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
source('normalDataParametricLikelihood.R')

nEquals100[2] = coverageProbability


### K=0.5n
breakLocs <- 0.5 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
source('normalDataParametricLikelihood.R')

nEquals100[3] = coverageProbability


### K=0.75n
breakLocs <- 0.75 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
source('normalDataParametricLikelihood.R')

nEquals100[4] = coverageProbability



################################################################################
## n = 200
simLen <- 200
### K=0.2n
breakLocs <- 0.2 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
source('normalDataParametricLikelihood.R')

nEquals200[1] = coverageProbability


### K=0.3n
breakLocs <- 0.3 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
source('normalDataParametricLikelihood.R')

nEquals200[2] = coverageProbability


### K=0.5n
breakLocs <- 0.5 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
source('normalDataParametricLikelihood.R')

nEquals200[3] = coverageProbability


### K=0.75n
breakLocs <- 0.75 * simLen
simMatrix <- indepDatagen(simLen=simLen, numSims=numSims, breakLocations=breakLocs, param1=param1, param2=param2)
source('normalDataParametricLikelihood.R')

nEquals200[4] = coverageProbability



################################################################################
# Constructing the output table
results.data <- data.frame(
  breakLocs = c("K=0.2n", "K=0.3n", "K=0.5n", "K=0.75n"),
  nEquals50 = nEquals50,
  nEquals100 = nEquals100,
  nEquals200 = nEquals200
)

print(results.data)