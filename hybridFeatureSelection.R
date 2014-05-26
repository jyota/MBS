hybridFeatureSelection <- function(x, y, start="random",stopP, stopT2){
  # Implements stepwise hybrid feature selection
  # adapted from Darius Dziuda's 2010 "Data Mining for Genomics & Proteomics" book
  # x is dataframe or matrix of predictors/variables that vary per class
  # y is class vector or dataframe
  # start is choice of starting variable. Can be "random", "Hotelling", or "randomForest"
  # stopP will stop search at a specified number of variables
  # stopT2 will stop search at a certain T2 threshold 
  source("obtainBestInitial.R")
  source("forwardSelection.R")
  source("backwardOptimize.R")
  if(is.null(stopP) | is.null(stopT2)){
    cat("Must enter both stop criterion.\n")
    return(NULL)
  }
  pool = x
  #cat("Beginning search for ", stopP, " variables.\n")
  currentSet = NULL
  poolSize   = NCOL(x)
  markerSize = 0
  currentT2  = 0.0
  i          = 1
  while(markerSize < stopP && currentT2 < stopT2 && (markerSize <= (NROW(x)-NROW(unique(y))-2) | (NROW(x)-NROW(unique(y))-2) == 0) && markerSize < NCOL(x)){
    maxGainStep = 0.0
    markerSize = markerSize + 1
    if(markerSize == 1){
      init = obtainBestInitial(as.matrix(x),as.matrix(y),method=start)
      selectedVar = init$maxVar
      maxGainStep = init$maxGain    
    } else {
      stepForward = forwardSelection(as.matrix(currentSet),y,as.matrix(pool))
      selectedVar = stepForward$selectedVar
      maxGainStep = stepForward$maxGain
    }
    currentT2  = currentT2 + maxGainStep
    if(!is.null(currentSet)){
      currentSet = cbind(currentSet, subset(pool,select=c(selectedVar)))
    }else{
      currentSet = subset(pool,select=c(selectedVar))
    }
    pool = subset(pool,select=-c(selectedVar))
    poolSize = poolSize - 1
    if(markerSize > 2){
      stepBack = backwardOptimize(currentSet,y)
      minLossStep = stepBack$minLoss
      if(minLossStep < maxGainStep - 0.0001 && !is.null(stepBack$dropVar)){        
                       # adjustment above to maxGainStep is to avoid getting stuck on a variable with very slightly lowered minLossStep
	pool = cbind(pool, subset(currentSet,select=c(stepBack$dropVar)))
	currentSet = subset(currentSet,select=-c(stepBack$dropVar))
	poolSize = poolSize + 1
	markerSize = markerSize - 1
	currentT2 = currentT2 - minLossStep
      }
    }
    #cat("Iteration: ", i, ", Current stat calc: ", currentT2, "\nVars: ",colnames(currentSet), "\n")    
    i = i + 1
  }
  return(as.data.frame(currentSet))
 }
