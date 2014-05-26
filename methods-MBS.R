require(psych)

newMBS <- function(dataMatrix, classes, stopP = 5, stopT2 = 1000.0, reps = 1, initialSelection = "random", priors = NULL, proportionInBag = 0.632)
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
	mbs <- .MBS(dataMatrix = dataMatrix, classes = classes, stopP = stopP, stopT2 = stopT2, reps = reps, initialSelection = "random", priors = priors, proportionInBag = 0.632)
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

