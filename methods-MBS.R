require(psych)
require(MASS)

newMBS <- function(dataMatrix, classes, stopP = 5, stopT2 = 1000.0, reps = 1, initialSelection = "random", priors = NULL, proportionInBag = 0.632, searchWithReplacement = TRUE)
{
	dataMatrix = as.matrix(dataMatrix)
	classes = as.numeric(classes)
	if(is.null(priors)) priors <- rep(1.0 / length(unique(classes)), length(unique(classes)))
	if(length(classes) != nrow(dataMatrix)){
	  stop("Ensure length(classes) == nrow(dataMatrix) -- you may need to transpose your predictor matrix.\n")
	}
	if(length(priors) != length(unique(classes))){
	  stop("Ensure length(priors) == length(unique(classes)).\n")
	}
	mbs <- .MBS(dataMatrix = dataMatrix, classes = classes, stopP = stopP, stopT2 = stopT2, reps = reps, initialSelection = "random", priors = priors, proportionInBag = 0.632, searchWithReplacement = searchWithReplacement)
	return(mbs)
}

setMethod("getMBSIterationResults", signature(object = "MBS"),
	  function(object){
		object@iterationResults
	  }
)

setMethod("getMBSMetricResults", signature(object = "MBS"),
	  function(object){
		  list(avgAccuracy = object@avgAccuracy, avgSensitivity = object@avgSensitivity, avgSpecificity = object@avgSpecificity, avgT2 = object@avgT2)
	  }
)

setMethod("mbsMvar", signature(object = "MBS", selectedCols = "numeric", selectedRows = "numeric"),
	function(object, selectedCols, selectedRows){	
	# Implements multivariate regression -- produces test statistic that is also available with R's MANOVA, but we can calculate here for 1 explanatory variable unlike MANOVA.
	Y_ = matrix(ncol=NCOL(as.matrix(object@classes[selectedRows]))+1,nrow=NROW(as.matrix(object@classes[selectedRows])))
	Y_[ ,2] = as.matrix(object@classes[selectedRows])
	Y_[ ,1] = 1
	if(rcond(crossprod(Y_)) < .Machine$double.eps){
	  return(list(HotellingLawleyTrace=0.0))
	}else{
	  BETA = solve(crossprod(Y_)) %*% t(Y_) %*% object@dataMatrix[selectedRows, selectedCols]
	}
	
	X_ = Y_ %*% BETA
	ERROR_ = object@dataMatrix[selectedRows, selectedCols] - X_

	MEAN_X=matrix(ncol=NCOL(object@dataMatrix[selectedRows, selectedCols]),nrow=NROW(object@dataMatrix[selectedRows, selectedCols]))
	if (NCOL(object@dataMatrix[selectedRows, selectedCols])==1){
		MEAN_X[,1]=mean(object@dataMatrix[selectedRows, selectedCols])
	}else {
	for(k in 1:NCOL(object@dataMatrix[selectedRows, selectedCols])){
		MEAN_X[,k]=mean(object@dataMatrix[selectedRows, selectedCols][,k])
	}
	}

	SSCP_regression = crossprod(X_) - crossprod(MEAN_X)
	SSCP_residual   = crossprod(ERROR_)
	SSCP_total      = crossprod(object@dataMatrix[selectedRows, selectedCols]) - crossprod(MEAN_X)

	if(rcond(SSCP_residual) < .Machine$double.eps){
	  return(list(HotellingLawleyTrace=0.0))
	 }else{
	T2 = tr(SSCP_regression %*% solve(SSCP_residual))
    
	return(list(HotellingLawleyTrace=T2))
	  }
})

setMethod("mbsForwardSelection", signature(object = "MBS", selectedCols = "numeric", selectedRows = "numeric"),
	  function(object, selectedCols, selectedRows){
          # Forward selection maximizing Hotelling-Lawley trace
          # selectedCols is currently selected variables. 
	  # selectedRows is in-bag rows selected for modified bagging schema.
	  maxGain = 0.0  
	  currentT2 = mbsMvar(object, selectedCols, selectedRows)$HotellingLawleyTrace
 	pool = seq_len(NCOL(object@dataMatrix))
	pool = pool[-which(pool %in% selectedCols)]
  	for(i in seq_len(length(pool))){
    	  deltaT2 = mbsMvar(object, c(selectedCols, pool[i]), selectedRows)$HotellingLawleyTrace - currentT2
    
    	if(deltaT2 >= maxGain){
      	  maxGain = deltaT2
      	  selectedVar = pool[i]
    	}
  }  
  return(list(selectedVar=selectedVar,maxGain=maxGain))
})

setMethod("mbsBackwardOptimize", signature(object = "MBS", selectedCols = "numeric", selectedRows = "numeric"),
	  function(object, selectedCols, selectedRows){
          # Backward optimization with Hotelling-Lawley trace
	  # (this function returns a possible variable to remove and minLoss value--
	  #  if minLoss is less than maxGain from last forwardSelection variable may be dropped) 
	  dropVar=NULL
	  if(length(selectedCols)>2){
	      currentT2 = mbsMvar(object, selectedCols, selectedRows)$HotellingLawleyTrace      
	      minLoss = currentT2
	  for(i in seq_len(length(selectedCols))){
	      deltaT2 = currentT2 - mbsMvar(object, selectedCols[-i], selectedRows)$HotellingLawleyTrace
	    if(deltaT2 <= minLoss){
	      minLoss = deltaT2
	      dropVar = selectedCols[i]
	    }
	  }  
    	return(list(dropVar=dropVar,minLoss=minLoss))
    	} else{
      return(list(dropVar=NULL,minLoss=NULL))
      }
})

setMethod("mbsObtainBestInitial", signature(object = "MBS", selectedRows = "numeric"),
	function(object, selectedRows){
	# Obtain best initial variable for feature selection process.
	maxGain = 0.0
	maxVar    = NULL
	if(object@initialSelection == "Hotelling"){
		for(i in NCOL(object@dataMatrix)){
			currentTest = mbsMvar(object, selectedCols = i, selectedRows = selectedRows)
			if(currentTest$HotellingLawleyTrace > maxGain){			
				maxGain   = currentTest$HotellingLawleyTrace
				maxVar    = i
			}
		}
	}
	if(object@initialSelection == "random"){
		# Just pick a random column (returns Hotelling-Lawley trace stat)
		maxVar = sample(seq_len(NCOL(object@dataMatrix)),1)
		currentTest = mbsMvar(object, selectedCols = maxVar, selectedRows = selectedRows)
	        maxGain = currentTest$HotellingLawleyTrace
	}
	return(list(maxVar=maxVar,maxGain=maxGain))
})

setMethod("mbsHybridFeatureSelection", signature(object = "MBS", selectedRows = "numeric"),
	  function(object, selectedRows){
     	  if(is.null(object@stopP) | is.null(object@stopT2)){
    		stop("Must enter both stop criterion.\n")
    	  }
	pool = seq_len(NCOL(object@dataMatrix))
  	currentSet = NULL
  	poolSize   = length(pool)
  	markerSize = 0
  	currentT2  = 0.0
  	i          = 1
  	while(markerSize < object@stopP && currentT2 < object@stopT2 && (markerSize <= (length(selectedRows) - length(unique(object@classes)) - 2) | (length(selectedRows) - length(unique(object@classes)) - 2) == 0) && poolSize > 0){
    maxGainStep = 0.0
    markerSize = markerSize + 1
    if(markerSize == 1){
      init = mbsObtainBestInitial(object = object, selectedRows = selectedRows)
      selectedVar = init$maxVar
      maxGainStep = init$maxGain
    } else {
      stepForward = mbsForwardSelection(object = object, selectedCols = currentSet, selectedRows = selectedRows)
      selectedVar = stepForward$selectedVar
      maxGainStep = stepForward$maxGain
    }
    currentT2  = currentT2 + maxGainStep
    if(!is.null(currentSet)){
      currentSet = c(currentSet, selectedVar)
    }else{
      currentSet = selectedVar
    }
    pool = pool[-which(pool %in% selectedVar)]
    poolSize = poolSize - 1
    if(markerSize > 2){
      stepBack = mbsBackwardOptimize(object = object, selectedCols = currentSet, selectedRows = selectedRows)
      minLossStep = stepBack$minLoss
      if(minLossStep < maxGainStep - 0.0001 && !is.null(stepBack$dropVar)){        
        # adjustment above to maxGainStep is to avoid getting stuck on a variable with very slightly lowered minLossStep
	pool = c(pool, stepBack$dropVar)
	currentSet = currentSet[-which(currentSet %in% stepBack$dropVar)]
	poolSize = poolSize + 1
	markerSize = markerSize - 1
	currentT2 = currentT2 - minLossStep
      }
    }
    i = i + 1
  }
  return(currentSet)
})

setMethod("mbsRun", signature(object = "MBS"),
	  function(object){
		  # Implements modified bagging schema to obtain estimates related to feature selection
# x is data frame of independent variables
# y is data frame dependent variables
# rep is number of random samples/LDA classifiers to obtain estimates with
# start, stopP, stopT2 are used for hybridFeatureSelection procedure run prior to checking each classifier on Out of Bag samples
# proportion is proportion of x to be in bag, rest will be out of bag
		  cat("Running modified bagging schema for ", object@reps, " iterations...\n")
		  returnMatrix = matrix(nrow=object@reps,ncol=(2 + length(unique(object@classes))+object@stopP))
  		  returnMatrix = as.data.frame(returnMatrix)
		  colnames(returnMatrix) <- c('Index', 'T2', paste('Correct class ', unique(object@classes), sep = ""), paste('V',seq_len(object@stopP),sep=""))  		  
		  baggingProgressBar <- txtProgressBar(style=3)
  for(j in 1:rep){
      classFrame <- data.frame(id_seq = seq_len(length(object@classes)), class = object@classes)
      classFrame$selected <- FALSE
      for(z in unique(object@classes)){
	  inBagZ <- sample(classFrame[classFrame$class == z, ]$id_seq, size = round(nrow(classFrame[classFrame$class == z, ])*object@proportionInBag, 0), replace = FALSE)
      	  classFrame[inBagZ, ]$selected <- TRUE
      }
      inBagJ1 = data.frame(x,y=y, check.names=FALSE)
      rownames(inBagJ1) = rownames(x)
      inBagJ1 = inBagJ1[inBagJ1$y==0,]      
      inBagJ2 = data.frame(x,y=y, check.names=FALSE)
      rownames(inBagJ2) = rownames(x)
      inBagJ2 = inBagJ2[inBagJ2$y==1,]
      inBagJ1 = inBagJ1[sample(nrow(inBagJ1),size=round(nrow(inBagJ1)*proportion,0),replace=FALSE),]
      inBagJ2 = inBagJ2[sample(nrow(inBagJ2),size=round(nrow(inBagJ2)*proportion,0),replace=FALSE),]      
      fullInBag = rbind(inBagJ1,inBagJ2)
      # Need to ensure there will not be a inverse singularity below
      fullInBag = fullInBag[,which(abs(round(colSums(fullInBag),0))!=0)]
      OOBJ1   = data.frame(x,y=y, check.names=FALSE)
      rownames(OOBJ1) = rownames(x)
      OOBJ1   = OOBJ1[OOBJ1$y==0,]
      OOBJ2   = data.frame(x,y=y, check.names=FALSE)
      rownames(OOBJ2) = rownames(x)
      OOBJ2   = OOBJ2[OOBJ2$y==1,]
      OOBJ1   = subset(OOBJ1, !(rownames(OOBJ1) %in% rownames(inBagJ1)))
      OOBJ2   = subset(OOBJ2, !(rownames(OOBJ2) %in% rownames(inBagJ2)))
      fullOOB = rbind(OOBJ1, OOBJ2)
      #fullOOB = fullOOB[,which(round(colSums(fullOOB),0)!=0]
      #cat("Beginning feature selection-- run #: ", j, " in bag samples: ", nrow(fullInBag), " OOB samples: ", nrow(fullOOB), "\n")
      tmpDat = hybridFeatureSelection(as.matrix(fullInBag[,1:(NCOL(fullInBag)-1)]),as.matrix(fullInBag[,NCOL(fullInBag)]),start,stopP,stopT2) 
      #cat("Beginning LDA fit-- run #: ", j, " ")
      if(!is.null(priors)){
      tmpFit = lda(classes ~ .,data=data.frame(tmpDat,classes=fullInBag[,NCOL(fullInBag)],check.names=FALSE),prior=priors)
      }else{
      tmpFit = lda(classes ~ .,data=data.frame(tmpDat,classes=fullInBag[,NCOL(fullInBag)],check.names=FALSE))
      }
      q = data.frame(y=as.factor(fullOOB$y),predict=predict(tmpFit,fullOOB)$class)
      #cat(" showing ", NROW(q[q[,1]==1,]), " for class 2, ", NROW(q[q[,1]==0,]), " for class 1, ", NROW(q[q[,1]==q[,2] & q[,1]==0,]), " correct class 1,", NROW(q[q[,1]==q[,2] & q[,1]==1,]), " correct class 2,", NROW(q[q[,1]==q[,2] & q[,1]==0,])+NROW(q[q[,1]==q[,2] & q[,1]==1,]), " correctly classified OOB samples.\n")
      #print(q)
      repStats[j,1] = (NROW(q[q[,1]==q[,2] & q[,1]==1,])+NROW(q[q[,1]==q[,2] & q[,1]==0,]))/NROW(fullOOB)
      repStats[j,2] = NROW(q[q[,1]==q[,2] & q[,1]==1,])/NROW(q[q[,1]==1,])
      repStats[j,3] = NROW(q[q[,1]==q[,2] & q[,1]==0,])/NROW(q[q[,1]==0,])
      repStats[j,4] = mvar(X=as.matrix(tmpDat),Y=as.matrix(fullInBag[,NCOL(fullInBag)]))$HotellingLawleyTrace
      #cat("result in accuracy: ", repStats[j,1], " sensitivity: ", repStats[j,2], " specificity: ", repStats[j,3],"\n")
      varsStats[varsStats$Variable %in% colnames(tmpDat),]$Times_Selected = varsStats[varsStats$Variable %in% colnames(tmpDat),]$Times_Selected + 1
      if(progressBar==T){
        setTxtProgressBar(baggingProgressBar,j/rep)
      }
      if(repStats[j,1]==1){
	varsStats[varsStats$Variable %in% colnames(tmpDat),]$Perfect_Selected = varsStats[varsStats$Variable %in% colnames(tmpDat),]$Perfect_Selected + 1
      }
  }
  varsStats$Perc_Selected = varsStats$Times_Selected / sum(varsStats$Times_Selected)
  if(progressBar==T){
     cat('\n')
  }
  return(list(varsStats=varsStats,repStats=repStats))

})

