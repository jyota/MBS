setGeneric("getMBSIterationResults", function(object) standardGeneric("getMBSIterationResults"))

setGeneric("getMBSMetricResults", function(object) standardGeneric("getMBSMetricResults"))

setGeneric("mbsBackwardOptimize", function(object, selectedCols, selectedRows) standardGeneric("mbsBackwardOptimize"))

setGeneric("mbsObtainBestInitial", function(object, selectedRows) standardGeneric("mbsObtainBestInitial"))

setGeneric("mbsHybridFeatureSelection", function(object, selectedRows) standardGeneric("mbsHybridFeatureSelection"))

setGeneric("mbsRun", function(object) standardGeneric("mbsRun"))

setGeneric("MBS", function(dataMatrix, classes, ...) standardGeneric("MBS"))
