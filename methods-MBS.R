require(psych)
require(MASS)

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
	Y_ = matrix(ncol = NCOL(as.matrix(object@classes[selectedRows])) + 1, nrow = NROW(as.matrix(object@classes[selectedRows])))
	Y_[ ,2] = as.matrix(object@classes[selectedRows])
	Y_[ ,1] = 1
	if(rcond(crossprod(Y_)) < .Machine$double.eps){
	  return(list(HotellingLawleyTrace=0.0))
	}else{
	  BETA = solve(crossprod(Y_)) %*% t(Y_) %*% object@dataMatrix[selectedRows, selectedCols]
	}
	
	X_ = Y_ %*% BETA
	ERROR_ = object@dataMatrix[selectedRows, selectedCols] - X_

	MEAN_X = matrix(ncol = NCOL(object@dataMatrix[selectedRows, selectedCols]), nrow = NROW(object@dataMatrix[selectedRows, selectedCols]))
	if (NCOL(object@dataMatrix[selectedRows, selectedCols])==1){
		MEAN_X[,1] = mean(object@dataMatrix[selectedRows, selectedCols])
	}else {
	for(k in 1:NCOL(object@dataMatrix[selectedRows, selectedCols])){
		MEAN_X[,k] = mean(object@dataMatrix[selectedRows, selectedCols][,k])
	}
	}

	SSCP_regression = crossprod(X_) - crossprod(MEAN_X)
	SSCP_residual   = crossprod(ERROR_)
	SSCP_total      = crossprod(object@dataMatrix[selectedRows, selectedCols]) - crossprod(MEAN_X)

	if(rcond(SSCP_residual) < .Machine$double.eps){
	  return(list(HotellingLawleyTrace = 0.0))
	 }else{
	T2 = tr(SSCP_regression %*% solve(SSCP_residual))
    
	return(list(HotellingLawleyTrace = T2))
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
    return(list(selectedVar = selectedVar, maxGain = maxGain))
})

setMethod("mbsBackwardOptimize", signature(object = "MBS", selectedCols = "numeric", selectedRows = "numeric"),
	  function(object, selectedCols, selectedRows){
          # Backward optimization with Hotelling-Lawley trace
	  # (this function returns a possible variable to remove and minLoss value--
	  #  if minLoss is less than maxGain from last forwardSelection variable may be dropped) 
	  dropVar = NULL
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
    	  return(list(dropVar = dropVar, minLoss = minLoss))
    	} else {
          return(list(dropVar = NULL, minLoss = NULL))
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
	return(list(maxVar = maxVar, maxGain = maxGain))
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
  	  while(markerSize < object@stopP && currentT2 < object@stopT2 && (markerSize <= (length(selectedRows) - length(unique(object@classes)) - 2) | (length(selectedRows) - length(unique(object@classes)) - 2) == 0) && markerSize < ncol(object@dataMatrix)){
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
		cat("Running modified bagging schema for ", object@reps, " iterations...\n")
		if(object@assessOutOfBag == TRUE){
	  		returnMatrix = matrix(nrow=object@reps,ncol=(2 + length(unique(object@classes))*2+object@stopP))
  			returnMatrix = as.data.frame(returnMatrix)
			colnames(returnMatrix) <- c('Index', 'T2', paste('Correct class ', unique(object@classes), sep = ""), paste('N in class ', unique(object@classes), sep = ""), paste('V',seq_len(object@stopP),sep=""))  		  
		} else {
	  		returnMatrix = matrix(nrow=object@reps,ncol=(2 + object@stopP))
  			returnMatrix = as.data.frame(returnMatrix)
			colnames(returnMatrix) <- c('Index', 'T2', paste('V',seq_len(object@stopP),sep=""))  		  
		}
		if(object@searchWithReplacement == FALSE){
			origDm <- object@dataMatrix
		}
		
		baggingProgressBar <- txtProgressBar(style=3)
  		 for(j in 1:object@reps){
		  if(object@assessOutOfBag == TRUE){	 
			  classFrame <- data.frame(id_seq = seq_len(length(object@classes)), class = object@classes)
			  classFrame$selected <- FALSE
	    		  for(z in unique(object@classes)){
		  		  inBagZ <- sample(classFrame[classFrame$class == z, ]$id_seq, size = round(nrow(classFrame[classFrame$class == z, ])*object@proportionInBag, 0), replace = FALSE)
		      		  classFrame[inBagZ, ]$selected <- TRUE
		 	 }
	         	tmpSelected = mbsHybridFeatureSelection(object = object, selectedRows = classFrame[classFrame$selected == TRUE, ]$id_seq) 
	          	if(length(tmpSelected)>1){
		      		fitDf <- data.frame(object@dataMatrix[classFrame[classFrame$selected == TRUE, ]$id_seq, tmpSelected],classes=classFrame[classFrame$selected == TRUE, ]$class, check.names = FALSE)
	          	} else {
				fitDf <- data.frame(variable=object@dataMatrix[classFrame[classFrame$selected == TRUE, ]$id_seq, tmpSelected],classes = classFrame[classFrame$selected == TRUE, ]$class, check.names = FALSE)
				colnames(fitDf) <- c(colnames(object@dataMatrix)[tmpSelected], 'classes')
	          	}
	       		tmpFit = lda(classes ~ ., data = fitDf, prior = object@priors)
	      	        if(length(tmpSelected)>1){
	 	     	 	predictDf <- data.frame(object@dataMatrix[classFrame[classFrame$selected == FALSE, ]$id_seq, tmpSelected],check.names=FALSE)
	       		} else {
				predictDf <- data.frame(variable=object@dataMatrix[classFrame[classFrame$selected == FALSE, ]$id_seq, tmpSelected])
	 			colnames(predictDf) <- colnames(object@dataMatrix)[tmpSelected]
	         	}
	         	q = data.frame(y = classFrame[classFrame$selected == FALSE, ]$class, predict = predict(tmpFit, predictDf)$class)
			returnMatrix[j, ]$Index <- j
			returnMatrix[j, ]$T2 <- mbsMvar(object, selectedCols = tmpSelected, selectedRows = classFrame[classFrame$selected == FALSE, ]$id_seq)$HotellingLawleyTrace
			tmpC <- unique(object@classes)
			for(v in 1:length(tmpC)){
			      returnMatrix[j, 2 + v] <- nrow(q[q$y == tmpC[v] & q$predict == tmpC[v], ])
			      returnMatrix[j, 2 + length(tmpC) + v] <- nrow(q[q$y == tmpC[v], ])
	      		 }
		         returnMatrix[j, (2 + length(tmpC)*2 + 1):ncol(returnMatrix)] <- colnames(object@dataMatrix)[tmpSelected]
		    }else{
		       #Not assessing out of bag importance, so just select variables and calculate T2 against whole data matrix.
	               tmpSelected = mbsHybridFeatureSelection(object = object, selectedRows = seq_len(nrow(object@dataMatrix))) 
    		       returnMatrix[j, ]$Index <- j
	               returnMatrix[j, ]$T2 <- mbsMvar(object, selectedCols = tmpSelected, selectedRows = seq_len(nrow(object@dataMatrix)))$HotellingLawleyTrace
	               returnMatrix[j, 3:ncol(returnMatrix)] <- colnames(object@dataMatrix)[tmpSelected]
		    }
		    object@usedVars <- unique(c(object@usedVars, tmpSelected))
		    if(object@searchWithReplacement == FALSE & length(object@usedVars) > 0){
	  		object@dataMatrix <- object@dataMatrix[, -object@usedVars]
	  	    }
		    setTxtProgressBar(baggingProgressBar,j/object@reps)
		}
	   object@iterationResults <- returnMatrix
	   object@avgT2 <- mean(returnMatrix$T2)
	   if(object@assessOutOfBag == TRUE){
		tmpC <- unique(object@classes)
	   	accCalc <- rowSums(returnMatrix[, 3:(2+length(tmpC))]) / rowSums(returnMatrix[, (3+length(tmpC)):(((3+length(tmpC)) - 1)+length(tmpC))])
	   	object@avgAccuracy <- mean(accCalc)
	   }
	   cat("\n")
	   if(object@searchWithReplacement == FALSE){
		object@dataMatrix <- origDm
	   }	   
	   return(object)
})

setMethod("MBS", signature(dataMatrix = "matrix", classes = "numeric"), 
   function(dataMatrix, classes, stopP = 5, stopT2 = 1000.0, reps = 1, initialSelection = "random", priors = NULL, proportionInBag = 0.632, searchWithReplacement = TRUE, assessOutOfBag = TRUE)
	{
	classes = as.numeric(classes)
	if(is.null(priors)) priors <- rep(1.0 / length(unique(classes)), length(unique(classes)))
	if(length(classes) != nrow(dataMatrix)){
	  stop("Ensure length(classes) == nrow(dataMatrix) -- you may need to transpose your predictor matrix.\n")
	}
	if(length(priors) != length(unique(classes))){
	  stop("Ensure length(priors) == length(unique(classes)).\n")
	}
	if(stopP > ncol(dataMatrix)){
	  stop("stopP > ncol(dataMatrix)! Aborting.\n")
	}
	mbs <- .MBS(dataMatrix = dataMatrix, classes = classes, stopP = stopP, stopT2 = stopT2, reps = reps, initialSelection = initialSelection, priors = priors, proportionInBag = proportionInBag, searchWithReplacement = searchWithReplacement, assessOutOfBag = assessOutOfBag)
	return(mbsRun(mbs))
})

setMethod("MBS", signature(dataMatrix = "data.frame", classes = "numeric"),
	  function(dataMatrix, classes, ...){ 
		  MBS(as.matrix(dataMatrix), classes = classes, ...) 
})

