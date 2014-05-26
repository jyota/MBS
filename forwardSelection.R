forwardSelection <- function(x, y, pool){
  # Forward selection maximizing Hotelling-Lawley trace
  # x is current dataframe or matrix of predictors or variables varying on class
  # y is class
  # pool is pool of available variables to add to x 
  # Adapted from Darius Dziuda's 2010 "Data Mining for Genomics & Proteomics" book
  source("multivariate-regression.R")

  maxGain = 0.0  
  currentT2 = mvar(as.matrix(x),y)$HotellingLawleyTrace
    
  for(i in 1:NCOL(pool)){
    deltaT2 = mvar(as.matrix(cbind(pool[,i],x)),y)$HotellingLawleyTrace - currentT2
    
    if(deltaT2 >= maxGain){
      maxGain = deltaT2
      selectedVar = i
    }
 }   
  return(list(selectedVar=selectedVar,maxGain=maxGain))
}
 
