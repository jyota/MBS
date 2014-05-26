modifiedBagging <- function(x, y, rep=1000, proportion=0.632, start="random",stopP,stopT2,priors=NULL)
{
# Implements modified bagging schema to obtain estimates related to feature selection
# x is data frame of independent variables
# y is data frame dependent variables
# rep is number of random samples/LDA classifiers to obtain estimates with
# start, stopP, stopT2 are used for hybridFeatureSelection procedure run prior to checking each classifier on Out of Bag samples
# proportion is proportion of x to be in bag, rest will be out of bag
  source("hybridFeatureSelection.R")
  require(MASS)
  repStats = as.data.frame(matrix(ncol=4,nrow=rep))
  colnames(repStats) = c("Accuracy","Sensitivity","Specificity", "T2")
  varsStats = as.data.frame(matrix(ncol=4, nrow=ncol(x)))
  colnames(varsStats) = c("Variable", "Times_Selected","Perc_Selected","Perfect_Selected")
  varsStats[,1] = as.matrix(colnames(x))
  varsStats[,2] = 0.0
  varsStats[,3] = 0.0
  varsStats[,4] = 0.0

  for(j in 1:rep){
      # for now, handles only two class problem
      inBagJ1 = data.frame(x,y=y, check.names=FALSE)
      rownames(inBagJ1) = rownames(x)
      inBagJ1 = inBagJ1[inBagJ1$y==0,]      
      inBagJ2 = data.frame(x,y=y, check.names=FALSE)
      rownames(inBagJ2) = rownames(x)
      inBagJ2 = inBagJ2[inBagJ2$y==1,]
      inBagJ1 = inBagJ1[sample(nrow(inBagJ1),size=round(nrow(inBagJ1)*proportion,0),replace=FALSE),]
      inBagJ2 = inBagJ2[sample(nrow(inBagJ2),size=round(nrow(inBagJ2)*proportion,0),replace=FALSE),]      
      fullInBag = rbind(inBagJ1,inBagJ2)
      # Need to ensure there will not be a inverse singularity below
      fullInBag = fullInBag[,which(abs(round(colSums(fullInBag),0))!=0)]
      OOBJ1   = data.frame(x,y=y, check.names=FALSE)
      rownames(OOBJ1) = rownames(x)
      OOBJ1   = OOBJ1[OOBJ1$y==0,]
      OOBJ2   = data.frame(x,y=y, check.names=FALSE)
      rownames(OOBJ2) = rownames(x)
      OOBJ2   = OOBJ2[OOBJ2$y==1,]
      OOBJ1   = subset(OOBJ1, !(rownames(OOBJ1) %in% rownames(inBagJ1)))
      OOBJ2   = subset(OOBJ2, !(rownames(OOBJ2) %in% rownames(inBagJ2)))
      fullOOB = rbind(OOBJ1, OOBJ2)
      #fullOOB = fullOOB[,which(round(colSums(fullOOB),0)!=0]
      cat("Beginning feature selection-- run #: ", j, " in bag samples: ", nrow(fullInBag), " OOB samples: ", nrow(fullOOB), "\n")
      tmpDat = hybridFeatureSelection(as.matrix(fullInBag[,1:(NCOL(fullInBag)-1)]),as.matrix(fullInBag[,NCOL(fullInBag)]),start,stopP,stopT2) 
      cat("Beginning LDA fit-- run #: ", j, " ")
      if(!is.null(priors)){
      tmpFit = lda(classes ~ .,data=data.frame(tmpDat,classes=fullInBag[,NCOL(fullInBag)],check.names=FALSE),prior=priors)
      }else{
      tmpFit = lda(classes ~ .,data=data.frame(tmpDat,classes=fullInBag[,NCOL(fullInBag)],check.names=FALSE))
      }
      q = data.frame(y=as.factor(fullOOB$y),predict=predict(tmpFit,fullOOB)$class)
      cat(" showing ", NROW(q[q[,1]==1,]), " for class 2, ", NROW(q[q[,1]==0,]), " for class 1, ", NROW(q[q[,1]==q[,2] & q[,1]==0,]), " correct class 1,", NROW(q[q[,1]==q[,2] & q[,1]==1,]), " correct class 2,", NROW(q[q[,1]==q[,2] & q[,1]==0,])+NROW(q[q[,1]==q[,2] & q[,1]==1,]), " correctly classified OOB samples.\n")
      print(q)
      repStats[j,1] = (NROW(q[q[,1]==q[,2] & q[,1]==1,])+NROW(q[q[,1]==q[,2] & q[,1]==0,]))/NROW(fullOOB)
      repStats[j,2] = NROW(q[q[,1]==q[,2] & q[,1]==1,])/NROW(q[q[,1]==1,])
      repStats[j,3] = NROW(q[q[,1]==q[,2] & q[,1]==0,])/NROW(q[q[,1]==0,])
      repStats[j,4] = mvar(X=as.matrix(tmpDat),Y=as.matrix(fullInBag[,NCOL(fullInBag)]))$HotellingLawleyTrace
      cat("result in accuracy: ", repStats[j,1], " sensitivity: ", repStats[j,2], " specificity: ", repStats[j,3],"\n")
      varsStats[varsStats$Variable %in% colnames(tmpDat),]$Times_Selected = varsStats[varsStats$Variable %in% colnames(tmpDat),]$Times_Selected + 1
      if(repStats[j,2] > 0.7 & repStats[j,3] > 0.7){
	varsStats[varsStats$Variable %in% colnames(tmpDat),]$Perfect_Selected = varsStats[varsStats$Variable %in% colnames(tmpDat),]$Perfect_Selected + 1
      }
  }
  varsStats$Perc_Selected = varsStats$Times_Selected / sum(varsStats$Times_Selected)
  return(list(varsStats=varsStats,repStats=repStats))
}
