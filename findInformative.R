# Generate a sequence of alternative biomarkers, where variables used in previous iteration
# are removed from pool of possible variables, for a given P variables. (entire training set)

findInformative <- function(x, y, rep=100, proportion=0.632, start="random",stopP,stopT2,priors=NULL){
  source("hybridFeatureSelection.R")
  require(MASS)
  
  returnMatrix = matrix(nrow=rep,ncol=(2+stopP))
   
  cat("Finding informative set of genes for ", rep, " iterations...\n")
  baggingProgressBar <- txtProgressBar(style=3)
  result = hybridFeatureSelection(as.matrix(x),as.matrix(y),start=start,stopP=stopP,stopT2=stopT2)
  finals = colnames(result)
  returnMatrix[1,] = c(1,mvar(X=as.matrix(result),Y=as.matrix(classes))$HotellingLawleyTrace, finals)
  for(i in 2:rep){    
   #cat("Iteration ", i, "..., ", length(finals), " variables removed, ", NCOL(x[,-which(colnames(x) %in% finals)]), " remaining.\n")
   if(NCOL(as.matrix(x[,-which(colnames(x) %in% finals)]))>stopP){
     result = hybridFeatureSelection(as.matrix(x[,-which(colnames(x) %in% finals)]),as.matrix(y),start=start,stopP=stopP,stopT2=stopT2)
     finals = c(colnames(result),finals)
     returnMatrix[i,] = c(i,mvar(X=as.matrix(result),Y=as.matrix(classes))$HotellingLawleyTrace, colnames(result))
     setTxtProgressBar(baggingProgressBar,i/rep)
   }else{
    i <- rep
   }
   #cat(returnMatrix[i,2], " is T2 for this iteration.\n")
  }
  returnMatrix = as.data.frame(returnMatrix)
  colnames(returnMatrix) <- c('Index','T2',paste('V',seq(1:stopP),sep=""))
  cat('\n')
  return(returnMatrix)
}
