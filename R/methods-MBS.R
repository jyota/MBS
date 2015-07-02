setMethod("summary", signature(object = "MBS"),
	  function(object, variableUsedAccuracy = 1.0, numberOfAccurateVars = 10){
		if(object@assessOutOfBag == TRUE){
			tmpC <- unique(object@classes)
	   		accCalc <- rowSums(object@iterationResults[, 3:(2+length(tmpC))]) / rowSums(object@iterationResults[, (3+length(tmpC)):(((3+length(tmpC)) - 1)+length(tmpC))])
	   		perfectVars <- object@iterationResults[accCalc >= variableUsedAccuracy, which(grepl("V", colnames(object@iterationResults)))]
			showVars <- table(unlist(perfectVars))
			showVars <- showVars[order(showVars, decreasing = TRUE)]
			cat("Modified Bagging Schema ensemble summary\n",
		    		"Iterations: ", object@reps, "\n",
		    		"Average T2 statistic: ", object@avgT2, "\n",
		    		"Average accuracy: ", object@avgAccuracy, "\n",
		    		"Table of top ", numberOfAccurateVars, " used variables selected into classifiers with ", variableUsedAccuracy, " accuracy: \n") 
		    		print(showVars[1:numberOfAccurateVars])
		} else {
			cat("Modified Bagging Schema ensemble summary\n",
		    		"Iterations: ", object@reps, "\n",
		    		"Average T2 statistic: ", object@avgT2, "\n")

		}
	  }
)

setMethod("show", signature(object = "MBS"),
	 function(object){
		cat("Results of Modified Bagging Schema iterations:\n")
	 	object@iterationResults	
	 }
)

setMethod("mbsObtainBestInitial", signature(object = "MBS", selectedRows = "numeric"),
	function(object, selectedRows){
	# Obtain best initial variable for feature selection process.
	maxGain = 0.0
	maxVar    = NULL
	if(object@initialSelection == "Hotelling"){
		for(i in NCOL(object@dataMatrix)){
			if(length(unique(object@dataMatrix[selectedRows, i])) > 1){
				currentTest = mbsMvarR(as.matrix(object@dataMatrix[selectedRows, i]), object@classes[selectedRows])
				if(currentTest > maxGain){			
					maxGain   = currentTest
					maxVar    = i
				}
			}
		}
	}
	if(object@initialSelection == "random"){
		# Just pick a random column (returns Hotelling-Lawley trace stat)
		while(1 == 1){
		      maxVar <- sample(seq_len(NCOL(object@dataMatrix)),1)
		      if(length(unique(object@dataMatrix[selectedRows, maxVar])) > 1){
			maxGain <- mbsMvarR(as.matrix(object@dataMatrix[selectedRows, maxVar]), object@classes[selectedRows])
			break
		      }
		}
	}
	return(list(maxVar = maxVar, maxGain = maxGain))
})

setMethod("mbsHybridFeatureSelection", signature(object = "MBS", selectedRows = "numeric"),
	 function(object, selectedRows){
     	  if(is.null(object@stopP) | is.null(object@stopT2) | is.null(selectedRows)){
    		stop("Must enter both stop criterion and selected rows.\n")
    	  }
	  pool = seq_len(NCOL(object@dataMatrix))
  	  currentSet = NULL
  	  poolSize   = length(pool)
  	  markerSize = 0
  	  currentT2  = 0.0
  	  i          = 1
	  minLossStep = NA
	  maxGainStep = NA
  	  while(markerSize < object@stopP && currentT2 < object@stopT2 && (markerSize <= (length(selectedRows) - length(unique(object@classes)) - 2) | (length(selectedRows) - length(unique(object@classes)) - 2) == 0) && markerSize < ncol(object@dataMatrix)){
	  maxGainStep = 0.0
    	  markerSize = markerSize + 1
    	  if(markerSize == 1){
      		init = mbsObtainBestInitial(object = object, selectedRows = selectedRows)
      		selectedVar = init$maxVar
      		maxGainStep = init$maxGain
    	  } else {
      		stepForward = mbsForwardSelection(object@dataMatrix, selectedCols = currentSet - 1, selectedRows = selectedRows - 1, classes = object@classes)
		if(stepForward[1] != -1){
       	        	selectedVar = stepForward[1] + 1
      			maxGainStep = stepForward[2]
		}
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
      		stepBack = mbsBackwardOptimize(x = object@dataMatrix, classes = object@classes, selectedCols = currentSet - 1, selectedRows = selectedRows - 1)
	  	if(stepBack[1] != -1){
		        minLossStep = stepBack[2]
			dropVar = stepBack[1] +	1
		}
	       if(minLossStep < maxGainStep - 0.0001 && dropVar > 0){        
	       # adjustment above to maxGainStep is to avoid getting stuck on a variable with very slightly lowered minLossStep
			pool = c(pool, dropVar)
			currentSet = currentSet[-which(currentSet %in% dropVar)]
			poolSize = poolSize + 1
			markerSize = markerSize - 1
			currentT2 = currentT2 - minLossStep
	      }
	  }
	  i = i + 1
	  }

	return(currentSet)
})

setMethod("mbsRun", signature(object = "MBS", showProgress = "logical"),
	  function(object, showProgress){
		# Implements modified bagging schema to obtain estimates related to feature selection
		if(showProgress == TRUE) cat("Running modified bagging schema for ", object@reps, " iterations...\n")
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
		if(showProgress == TRUE) baggingProgressBar <- txtProgressBar(style=3)
  		 for(j in 1:object@reps){
		  if(object@assessOutOfBag == TRUE){	 
			  classFrame <- data.frame(id_seq = seq_len(length(object@classes)), class = object@classes)
			  classFrame$selected <- FALSE
	    		  for(z in unique(object@classes)){
		  		  inBagZ <- sample(classFrame[classFrame$class == z, ]$id_seq, size = round(nrow(classFrame[classFrame$class == z, ])*object@proportionInBag, 0), replace = FALSE)
		      		  classFrame[inBagZ, ]$selected <- TRUE
		 	 }
	         	tmpRes <- capture.output(tmpSelected <- mbsHybridFeatureSelection(object = object, selectedRows = classFrame[classFrame$selected == TRUE, ]$id_seq))
  			  if(length(tmpSelected)>1){
		      		fitDf <- data.frame(object@dataMatrix[classFrame[classFrame$selected == TRUE, ]$id_seq, tmpSelected], classes = classFrame[classFrame$selected == TRUE, ]$class, check.names = FALSE)
	          	} else {
				fitDf <- data.frame(variable=object@dataMatrix[classFrame[classFrame$selected == TRUE, ]$id_seq, tmpSelected], classes = classFrame[classFrame$selected == TRUE, ]$class, check.names = FALSE)
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
			returnMatrix[j, ]$T2 <- mbsMvarR(object@dataMatrix[classFrame[classFrame$selected == TRUE, ]$id_seq, tmpSelected], object@classes[classFrame[classFrame$selected == TRUE, ]$id_seq])
			tmpC <- unique(object@classes)
			for(v in 1:length(tmpC)){
			      returnMatrix[j, 2 + v] <- nrow(q[q$y == tmpC[v] & q$predict == tmpC[v], ])
			      returnMatrix[j, 2 + length(tmpC) + v] <- nrow(q[q$y == tmpC[v], ])
	      		 }
		         returnMatrix[j, (2 + length(tmpC)*2 + 1):ncol(returnMatrix)] <- colnames(object@dataMatrix)[tmpSelected]
		    }else{
		       #Not assessing out of bag importance, so just select variables and calculate T2 against whole data matrix.
	               tmpRes <- capture.output(tmpSelected <- mbsHybridFeatureSelection(object = object, selectedRows = seq_len(nrow(object@dataMatrix)))) 
    		       returnMatrix[j, ]$Index <- j
	               returnMatrix[j, ]$T2 <- mbsMvarR(object@dataMatrix[, tmpSelected], object@classes)
	               returnMatrix[j, 3:ncol(returnMatrix)] <- colnames(object@dataMatrix)[tmpSelected]
		    }
		    object@usedVars <- unique(c(object@usedVars, tmpSelected))
		    if(object@searchWithReplacement == FALSE & length(object@usedVars) > 0){
	  		object@dataMatrix <- data.matrix(object@dataMatrix[, -object@usedVars])
	  	    }
		    if(showProgress == TRUE) setTxtProgressBar(baggingProgressBar,j/object@reps)
		}
	   object@iterationResults <- returnMatrix
	   object@avgT2 <- mean(returnMatrix$T2)
	   if(object@assessOutOfBag == TRUE){
		tmpC <- unique(object@classes)
	   	accCalc <- rowSums(returnMatrix[, 3:(2+length(tmpC))]) / rowSums(returnMatrix[, (3+length(tmpC)):(((3+length(tmpC)) - 1)+length(tmpC))])
	   	object@avgAccuracy <- mean(accCalc)
	   }
	   if(showProgress == TRUE) cat("\n")
	   if(object@searchWithReplacement == FALSE){
		object@dataMatrix <- origDm
	   }	   
	   return(object)
})

setMethod("MBS", signature(dataMatrix = "matrix", classes = "numeric"), 
   function(dataMatrix, classes, stopP = NULL, stopT2 = NULL, reps = NULL, initialSelection = "random", priors = NULL, proportionInBag = 0.632, searchWithReplacement = TRUE, assessOutOfBag = TRUE, showProgress = TRUE)
	{
	classes = as.numeric(classes)
	if(is.null(stopP)){
	  warning("No stopP value assigned, setting to maximum number of data matrix columns minus 1.\n")
	  stopP <- ncol(dataMatrix) - 1
        }	  
	if(is.null(stopT2)){
	  warning("No stopT2 value assigned, setting to 1000.\n")
	  stopT2 <- 1000.0	
	}
	if(is.null(reps)){
	  warning("No number of reps assigned, setting to 1.\n")
	  reps <- 1	
	}	
	if(is.null(priors)) priors <- rep(1.0 / length(unique(classes)), length(unique(classes)))
	if(length(classes) != nrow(dataMatrix)){
	  stop("Ensure length(classes) == nrow(dataMatrix) -- you may need to transpose your predictor matrix.\n")
	}
	if(length(priors) != length(unique(classes))){
	  stop("Ensure length(priors) == length(unique(classes)).\n")
	}
	uniqueLength <- apply(dataMatrix, 2, function(x) length(unique(x)))
	dataMatrix <- subset(dataMatrix, select = uniqueLength > 1)
	if(stopP > ncol(dataMatrix)){
	  stop("stopP > ncol(dataMatrix)! Aborting.\n")
	}
	if(searchWithReplacement == FALSE & stopP*reps > ncol(dataMatrix)){
	  stop("Request to perform MBS iterations, excluding variables from previous runs made, but not enough variables to meet number of reps requested. Aborting.\n")
	}
	mbs <- .MBS(dataMatrix = dataMatrix, classes = classes, stopP = stopP, stopT2 = stopT2, reps = reps, initialSelection = initialSelection, priors = priors, proportionInBag = proportionInBag, searchWithReplacement = searchWithReplacement, assessOutOfBag = assessOutOfBag)
	return(mbsRun(mbs, showProgress))
})

setMethod("MBS", signature(dataMatrix = "data.frame", classes = "numeric"),
	  function(dataMatrix, classes, ...){ 
		  MBS(dataMatrix = data.matrix(dataMatrix), classes = classes, ...) 
})

# Work in progress below for multicore bagging procedure (not yet functional)
setMethod("MBSparallel", signature(dataMatrix = "matrix", classes = "numeric", cores = "numeric"),
	  function(dataMatrix, classes, cores = detectCores(), stopP = NULL, stopT2 = NULL, reps = NULL, initialSelection = "random", priors = NULL, proportionInBag = 0.632, searchWithReplacement = TRUE, assessOutOfBag = TRUE, showProgress = TRUE){
		if(searchWithReplacement == FALSE){
			stop("searchWithReplace set to False is intended for a sequential search-- use MBS instead of MBSparallel.\n")
		}
  		classes = as.numeric(classes)
		if(is.null(stopP)){
		  warning("No stopP value assigned, setting to maximum number of data matrix columns minus 1.\n")
		  stopP <- ncol(dataMatrix) - 1
		}	  
		if(is.null(stopT2)){
		  warning("No stopT2 value assigned, setting to 1000.\n")
		  stopT2 <- 1000.0	
		}
		if(is.null(reps)){
		  warning("No number of reps assigned, setting to core.\n")
		  reps <- cores
		}	
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
		if((reps %% cores) > 0){
		  warning(paste("Number of reps / number of cores did not balance. Removing ", reps %% cores, " reps to balance.", sep = ""))
		}
		reps <- reps - (reps %% cores)
		reps = round(reps / cores, 0)

		if(searchWithReplacement == FALSE & stopP*reps > ncol(dataMatrix)){
		  stop("Request to perform MBS iterations, excluding variables from previous runs made, but not enough variables to meet number of reps requested. Aborting.\n")
		}
		if(showProgress == TRUE) cat(paste("Beginning ", reps, " iterations on each of ", cores, " cores (", reps*cores, " total)\n", sep = ""))
		mbs <- .MBS(dataMatrix = dataMatrix, classes = classes, stopP = stopP, stopT2 = stopT2, reps = reps, initialSelection = initialSelection, priors = priors, proportionInBag = proportionInBag, searchWithReplacement = searchWithReplacement, assessOutOfBag = assessOutOfBag)
		cl <- makeCluster(cores)
		clusterEvalQ(cl, {
			       library(MASS)
			       library(Rcpp)
			       library(RcppArmadillo)
			       source("Classes.R")
			       source("AllGenerics.R")
			       source("methods-MBS.R")
			       NULL })
		clusterExport(cl, "mbs", environment())
		mbsRes <- clusterEvalQ(cl, mbsRun(mbs, FALSE))
		stopCluster(cl)
		mbs@iterationResults <- mbsRes[[1]]@iterationResults
		mbs@usedVars <- mbsRes[[1]]@usedVars
		for(i in 2:length(mbsRes)){
			mbs@iterationResults <- rbind(mbs@iterationResults, mbsRes[[i]]@iterationResults)
			mbs@usedVars <- unique(c(mbs@usedVars, mbsRes[[i]]@usedVars))
		}
		mbs@avgT2 <- mean(mbs@iterationResults$T2)
		if(mbs@assessOutOfBag == TRUE){
			tmpC <- unique(mbs@classes)
			accCalc <- rowSums(mbs@iterationResults[, 3:(2+length(tmpC))]) / rowSums(mbs@iterationResults[, (3+length(tmpC)):(((3+length(tmpC)) - 1)+length(tmpC))])
			mbs@avgAccuracy <- mean(accCalc)
		}
		return(mbs)
})

