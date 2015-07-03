#' An S4 class for a modified bagging schema 
#' 
#' Instances should be created by calling the MBS method except for debug purposes.
#' @slot dataMatrix A matrix of data 
#' @slot classes Numeric label of class for each matrix row (class labels for training & assessing classifiers)
#' @slot stopP Number of variables stopping criterion
#' @slot stopT2 Hotelling-Lawley Trace statistic stopping criterion
#' @slot reps Bootstrap repetitions to use
#' @slot initialSelection Initial selection type (currently supports "random")
#' @slot proportionInBag Proportion of matrix rows as in-bag data for each bootstrap repetition
#' @slot priors Prior probability of classes to use for LDA classifiers
#' @slot avgAccuracy After mbsRun has completed, stores average accuracy of classifiers on Out-of-Bag data
#' @slot avgT2 After mbsRun has completed, stores average Hotelling-Lawley trace statistic value
#' @slot iterationResults After mbsRun has completed, has results for each bootstrap repetition
#' @slot assessOutOfBag Whether accuracy was assessed on Out-of-Bag samples
#' @slot searchWithReplacement After each bootstrap repetition, should full pool of columns be available for selection or should previously used columns be ignored
#' @slot usedVars Column numbers of variables selected as final features for any bootstrap repitition
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

