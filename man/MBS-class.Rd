% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/Classes.R
\docType{class}
\name{MBS-class}
\alias{.MBS}
\alias{MBS-class}
\title{An S4 class for a modified bagging schema}
\description{
Instances should be created by calling the MBS method except for debug purposes.
}
\section{Slots}{

\describe{
\item{\code{dataMatrix}}{A matrix of data}

\item{\code{classes}}{Numeric label of class for each matrix row (class labels for training & assessing classifiers)}

\item{\code{stopP}}{Number of variables stopping criterion}

\item{\code{stopT2}}{Hotelling-Lawley Trace statistic stopping criterion}

\item{\code{reps}}{Bootstrap repetitions to use}

\item{\code{initialSelection}}{Initial selection type (currently supports "random")}

\item{\code{proportionInBag}}{Proportion of matrix rows as in-bag data for each bootstrap repetition}

\item{\code{priors}}{Prior probability of classes to use for LDA classifiers}

\item{\code{avgAccuracy}}{After mbsRun has completed, stores average accuracy of classifiers on Out-of-Bag data}

\item{\code{avgT2}}{After mbsRun has completed, stores average Hotelling-Lawley trace statistic value}

\item{\code{iterationResults}}{After mbsRun has completed, has results for each bootstrap repetition}

\item{\code{assessOutOfBag}}{Whether accuracy was assessed on Out-of-Bag samples}

\item{\code{searchWithReplacement}}{After each bootstrap repetition, should full pool of columns be available for selection or should previously used columns be ignored}

\item{\code{usedVars}}{Column numbers of variables selected as final features for any bootstrap repitition}
}}

