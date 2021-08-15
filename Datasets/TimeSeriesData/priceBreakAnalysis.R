#### Souring the code to run the LR method on this, and also to run the MIC method on this data
source('LikelihoodRatio/normalDataParametricLikelihood.R')
source('SIC/normalIndepMICBreakEstimator.R')
source('LikelihoodRatio/normalDataConditionalLikelihood.R')

priceData <- read.csv(file='Datasets/TimeSeriesData/prices.csv')

cornMatrix <- matrix(priceData$Corn, nrow=1, ncol=length(priceData$Corn))
soybeanMatrix <- matrix(priceData$Soybean, nrow=1, ncol=length(priceData$Soybean))

# Corn results for the two methods
ILR_Corn_Results <- normalIndepLRCalc(cornMatrix, alpha=0.05)
MIC_Corn_Results <- normIndepMICCalc(cornMatrix, alpha=0.05)
CLR_Corn_Results <- normalCondLRCalc(cornMatrix, mu=mean(cornMatrix), sigma=sd(cornMatrix), alpha=0.05)

# Soybean results for the two methods
ILR_Soybean_Results <- normalIndepLRCalc(soybeanMatrix, alpha=0.05)
MIC_Soybean_Results <- normIndepMICCalc(soybeanMatrix, alpha=0.05)
CLR_Soybean_Results <- normalCondLRCalc(soybeanMatrix, mu=mean(soybeanMatrix), sigma=sd(soybeanMatrix), alpha=0.05)


