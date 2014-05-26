backwardOptimize <- function(x, y){
  # Backward optimization with Hotelling-Lawley trace
  # x is current dataframe or matrix of predictors or variables varying on class
  # (this function returns a possible variable to remove and minLoss value--
  #  if minLoss is less than maxGain from last forwardSelection variable may be dropped) 
  # y is class used)
  # Adapted from Darius Dziuda's 2010 "Data Mining for Genomics & Proteomics" book 
  source("multivariate-regression.R")
  dropVar=NULL
  if(NCOL(x)>2){
      currentT2 = mvar(as.matrix(x),y)$HotellingLawleyTrace      
    minLoss = currentT2
  for(i in 1:NCOL(x)){
      deltaT2 = currentT2 - mvar(as.matrix(x[,-i]),y)$HotellingLawleyTrace
    if(deltaT2 <= minLoss){
      minLoss = deltaT2
      dropVar = i
    }
  }  
  return(list(dropVar=dropVar,minLoss=minLoss))
  }
  else{
    cat("Predictor set is already 2 variables or less. Cancelling.\n")
    return(list(dropVar=NULL,minLoss=NULL))
  }
}
 
