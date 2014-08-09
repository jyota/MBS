\name{MBS-class}
\docType{class}
\alias{MBS-class}
\alias{show,MBS-method}
\alias{summary,MBS-method}
\title{Modified Bagging Schema - class}
\description{
An S4 class for storing results of iterations of modified bagging schema experiments.
}
\section{MBS Object Components}{
MBS objects are created as the result of a call to \code{MBS(dataMatrix = "matrix", classes = "numeric", ...)}.
This will create the object and run the modified bagging schema experiments for the chosen number of iterations.

Objects of the MBS class contain the following components:
  \describe{
    \item{\code{dataMatrix:}}{ matrix of the data used in the modified bagging experiments}
    \item{\code{classes:}}{ numeric vector of class labels }
    \item{\code{stopP:}}{ numeric stopping count of variables P to use for each iteration of feature selection in experiments }
    \item{\code{stopT2:}}{ numeric stopping value of T2 for each iteration of feature selection in experiments }
    \item{\code{reps:}}{ numeric of modified bagging experiment iterations performed/to perform }
    \item{\code{initialSelection:}}{ characters indicating type of initial selection in feature selection. can be "random" or "Hotelling" }
    \item{\code{proportionInBag:}}{ numeric proportion of records to keep in bag for train LDA models in modified bagging experiments }
    \item{\code{priors:}}{ numeric vector of priors assigned for classes in LDA models }
    \item{\code{avgAccuracy:}}{ numeric average out-of-bag accuracy observed in modified bagging schema experiments }
    \item{\code{avgT2:}}{ numeric average Lawley-Hotelling T2 observed for in-bag data during modified bagging schema experiments }
    \item{\code{iterationResults:}}{ data frame of results of each iteration of modified bagging schema experiments }
    \item{\code{assessOutOfBag:}}{ logical of whether to assess out of bag accuracy or just do feature selection on entire data set }
    \item{\code{searchWithReplacement:}}{ logical of whether to use variables from previous iterations in modified bagging schema experiments or not }
    \item{\code{usedVars:}}{ numeric vector tracking which columns have been used in modified bagging schema experiments if searchWithReplacement is set to false }
  }
}
\section{Methods}{
This class has methods that can be run directly for its objects if desired.
However, for end users these are typically unnecessary, as the methods MBS and MBSparallel can instead be used to initialize an MBS object and run the modified bagging schema experiments in one call.
  \describe{
    \item{getMBSIterationResults}{\code{signature(object = "MBS")}: ... }
    \item{getMBSMetricResults}{\code{signature(object = "MBS")}: ... }
    \item{mbsHybridFeatureSelection}{\code{signature(object = "MBS", selectedRows = "numeric")}: ... }
    \item{mbsObtainBestInitial}{\code{signature(object = "MBS", selectedRows = "numeric")}: ... }
    \item{mbsRun}{\code{signature(object = "MBS", showProgress = "logical")}: ... }
	 }
}
\references{
Dziuda DM (2010). \emph{Data mining for genomics and proteomics: Analysis of gene and protein expression data}. Hoboken, NJ: John Wiley & Sons, Inc.
}
\author{
Jyota Snyder <jyotas@my.ccsu.edu>. Modified bagging schema implemented was defined by Darius Dziuda.
}
\keyword{classes}

