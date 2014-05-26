require(methods)

# Modified bagging schema class
.MBS <- setClass("MBS", 
  representation(
	dataMatrix = "matrix",
	classes = "numeric",
	stopP = "numeric",
	stopT2 = "numeric",
	reps = "numeric",
	initialSelection = "character",
	proportionInBag = "numeric",
	priors = "numeric",
	avgAccuracy = "numeric",
	avgSensitivity = "numeric",
	avgSpecificity = "numeric",
	avgT2 = "numeric",
	iterationResults = "data.frame")
)

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

