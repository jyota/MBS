setGeneric("getMBSIterationResults", function(object) standardGeneric("getMBSIterationResults"))

setGeneric("getMBSMetricResults", function(object) standardGeneric("getMBSMetricResults"))

setGeneric("mbsObtainBestInitial", function(object, selectedRows) standardGeneric("mbsObtainBestInitial"))

setGeneric("mbsHybridFeatureSelection", function(object, selectedRows) standardGeneric("mbsHybridFeatureSelection"))

setGeneric("mbsRun", function(object, showProgress) standardGeneric("mbsRun"))

setGeneric("MBS", function(dataMatrix, classes, ...) standardGeneric("MBS"))

setGeneric("MBSparallel", function(dataMatrix, classes, cores, ...) standardGeneric("MBSparallel"))

