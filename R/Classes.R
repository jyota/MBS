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
	avgT2 = "numeric",
	iterationResults = "data.frame",
	assessOutOfBag = "logical",
	searchWithReplacement = "logical",
	usedVars = "numeric")
)

