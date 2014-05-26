mvar<-function(X, Y){	
	require(psych)
	# Implements multivariate regression -- produces two test statistics that are also available with R's MANOVA
	# X is matrix of variables that may vary per class (for example, gene expression levels)
	# Y is vector or matrix of class variable(s) (for example, a vector of binary classes specifying cancerous vs. non-cancerous)
	Y_ = matrix(ncol=NCOL(Y)+1,nrow=NROW(Y))
	Y_[,2:(NCOL(Y)+1)]=Y
	Y_[,1]=1
	if(rcond(crossprod(Y_)) < .Machine$double.eps){
	  return(list(HotellingLawleyTrace=0.0,WilksLambda=0.0))
	}else{
	BETA = solve(crossprod(Y_)) %*% t(Y_) %*% X
	}
	
	X_ = Y_ %*% BETA
	ERROR_ = X - X_

	MEAN_X=matrix(ncol=NCOL(X),nrow=NROW(X))
	if (NCOL(X)==1){
		MEAN_X[,1]=mean(X)
	}else {
	for(k in 1:NCOL(X)){
		MEAN_X[,k]=mean(X[,k])
	}
	}

	SSCP_regression = crossprod(X_) - crossprod(MEAN_X)
	SSCP_residual   = crossprod(ERROR_)
	SSCP_total      = crossprod(X) - crossprod(MEAN_X)

	WILKS = det(SSCP_residual)/det(SSCP_total)
	
	# W_CRITICAL = -((NCOL(Y)-NCOL(X)-1-((NCOL(Y)+1-NCOL(X))/2)))*log(WILKS)

	if(rcond(SSCP_residual) < .Machine$double.eps){
	  return(list(HotellingLawleyTrace=0.0,WilksLambda=0.0))
	 }else{
	T2 = tr(SSCP_regression %*% solve(SSCP_residual))
    
	return(list(HotellingLawleyTrace=T2,WilksLambda=WILKS))
	  }
}
