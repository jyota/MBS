\name{MBS-class}
\Rdversion{1.1}
\docType{class}
\alias{MBS-class}
\alias{getMBSIterationResults,MBS-method}
\alias{getMBSMetricResults,MBS-method}
\alias{mbsHybridFeatureSelection,MBS,numeric-method}
\alias{mbsObtainBestInitial,MBS,numeric-method}
\alias{mbsRun,MBS,logical-method}

\title{Class \code{"MBS"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MBS", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{dataMatrix}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{classes}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{stopP}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{stopT2}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{reps}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{initialSelection}:}{Object of class \code{"character"} ~~ }
    \item{\code{proportionInBag}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{priors}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{avgAccuracy}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{avgT2}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{iterationResults}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{assessOutOfBag}:}{Object of class \code{"logical"} ~~ }
    \item{\code{searchWithReplacement}:}{Object of class \code{"logical"} ~~ }
    \item{\code{usedVars}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{getMBSIterationResults}{\code{signature(object = "MBS")}: ... }
    \item{getMBSMetricResults}{\code{signature(object = "MBS")}: ... }
    \item{mbsHybridFeatureSelection}{\code{signature(object = "MBS", selectedRows = "numeric")}: ... }
    \item{mbsObtainBestInitial}{\code{signature(object = "MBS", selectedRows = "numeric")}: ... }
    \item{mbsRun}{\code{signature(object = "MBS", showProgress = "logical")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("MBS")
}
\keyword{classes}
