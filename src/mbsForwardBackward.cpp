// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector mbsMvarR(NumericMatrix x1, NumericVector classes1)
{
	try{
	arma::mat x = as<arma::mat>(x1);
	arma::vec classes = as<arma::vec>(classes1);
	arma::mat Y_ = arma::mat(classes.n_elem, 2);
	Y_.ones();
	Y_.col(1) = classes;

	arma::mat BETA = inv(trans(Y_) * Y_) * trans(Y_) * x;
	arma::mat X_ = Y_ * BETA;
	arma::mat ERROR_ = x - X_;
	arma::mat MEAN_X = arma::mat(x.n_rows, x.n_cols);
	for(int i = 0; i < x.n_cols; i++){
		MEAN_X.col(i) = arma::vec(x.n_rows).fill(arma::mean(arma::mean(x.col(i))));
	}

		arma::mat SSCP_regression = (trans(X_) * X_) - (trans(MEAN_X) * MEAN_X);
		arma::mat SSCP_residual = (trans(ERROR_) * ERROR_);
		arma::mat SSCP_total = (trans(x) * x) - (trans(MEAN_X) * MEAN_X);
	
		return wrap(trace(SSCP_regression * inv(SSCP_residual)));
	}catch(...){
		return(0.0);
	}
}

double mbsMvar(arma::mat x, arma::ivec classes)
{
	try{
	arma::mat Y_ = arma::mat(classes.n_elem, 2);
	Y_.ones();
	Y_.col(1) = arma::conv_to< arma::vec >::from(classes);

	arma::mat BETA = inv(trans(Y_) * Y_) * trans(Y_) * x;
	arma::mat X_ = Y_ * BETA;
	arma::mat ERROR_ = x - X_;
	arma::mat MEAN_X = arma::mat(x.n_rows, x.n_cols);
	for(int i = 0; i < x.n_cols; i++){
		MEAN_X.col(i) = arma::vec(x.n_rows).fill(arma::mean(arma::mean(x.col(i))));
	}

	arma::mat SSCP_regression = (trans(X_) * X_) - (trans(MEAN_X) * MEAN_X);
	arma::mat SSCP_residual = (trans(ERROR_) * ERROR_);
	arma::mat SSCP_total = (trans(x) * x) - (trans(MEAN_X) * MEAN_X);
	
	return trace(SSCP_regression * inv(SSCP_residual));	
	} catch(...){
		// If any kind of error, return 0.0 as T2 value.
		return 0.0;
	}
}

// [[Rcpp::export]]
NumericVector mbsForwardSelection(NumericMatrix x, NumericVector classes, NumericVector selectedCols, NumericVector selectedRows)
{
	double maxGain = 0.0;
	double deltaT2 = 0.0;
	int selectedVar = -1;
	arma::mat X = as<arma::mat>(x);
	arma::ivec Classes = as<arma::ivec>(classes);
	arma::ivec SelectedCols = as<arma::ivec>(selectedCols);
	arma::ivec SelectedRows = as<arma::ivec>(selectedRows);
	double currentT2 = mbsMvar(X.submat(arma::conv_to< arma::uvec >::from(SelectedRows), arma::conv_to< arma::uvec >::from(SelectedCols)), Classes(arma::conv_to< arma::uvec >::from(SelectedRows)));
	for(int i = 0; i < X.n_cols; i++){
		arma::uvec foundCol = arma::find(SelectedCols == i);
		if(foundCol.n_elem > 0){
			// Do nothing if already selected column selected.	
		} else {
			// Otherwise, combine this column with the others and get T2 value.
			arma::vec combCols = arma::vec(SelectedCols.n_elem + 1);
			combCols.fill(i);
			for(int j = 0; j < SelectedCols.n_elem; j++){
				combCols(j) = SelectedCols(j);
			}
			deltaT2 = mbsMvar(X.submat(arma::conv_to< arma::uvec >::from(SelectedRows), arma::conv_to< arma::uvec >::from(combCols)), Classes(arma::conv_to< arma::uvec >::from(SelectedRows))) - currentT2;

			if(deltaT2 >= maxGain & deltaT2 > 0.0){
				maxGain = deltaT2;
				selectedVar = i;
			}
		}
	}
	if(selectedVar != -1){
		arma::vec retCols = arma::vec(2);
		retCols(0) = selectedVar;
		retCols(1) = maxGain;
		return wrap(retCols);

	} else { 		
		arma::vec retCols = arma::vec(2);
		retCols(0) = -1;
		retCols(1) = -1;
		return(wrap(retCols));
	}
}


// [[Rcpp::export]]
NumericVector mbsBackwardOptimize(NumericMatrix x, NumericVector classes, NumericVector selectedCols, NumericVector selectedRows)
{
	int dropVar = -1;
	double minLoss = 0.0;
	double deltaT2 = 0.0;
	int tmpInc = 0;
	arma::mat X = as<arma::mat>(x);
	arma::ivec Classes = as<arma::ivec>(classes);
	arma::ivec SelectedCols = as<arma::ivec>(selectedCols);
	arma::ivec SelectedRows = as<arma::ivec>(selectedRows);
	if(SelectedCols.n_elem > 2){
		double currentT2 = mbsMvar(X.submat(arma::conv_to< arma::uvec >::from(SelectedRows), arma::conv_to< arma::uvec >::from(SelectedCols)), Classes(arma::conv_to< arma::uvec >::from(SelectedRows)));
		minLoss = currentT2;
		for(int i = 0; i < SelectedCols.n_elem; i++){
			arma::uvec reducedCols = arma::uvec(SelectedCols.n_elem - 1);	
			tmpInc = 0;
			for(int j = 0; j < SelectedCols.n_elem; j++){
				if(i != j){
					reducedCols(tmpInc) = SelectedCols(j);
					tmpInc++;
				}
			}
			deltaT2 = currentT2 - mbsMvar(X.submat(arma::conv_to< arma::uvec >::from(SelectedRows), reducedCols), Classes(arma::conv_to< arma::uvec >::from(SelectedRows)));
			if(deltaT2 <= minLoss){
				dropVar = SelectedCols(i);
				minLoss = deltaT2;
			}
		}
	}
	if(dropVar != -1){
		arma::vec retCols = arma::vec(2);
		retCols(0) = dropVar;
		retCols(1) = minLoss;
		return(wrap(retCols));
	} else { 		
		arma::vec retCols = arma::vec(2);
		retCols(0) = -1;
		retCols(1) = -1;
		return(wrap(retCols));
	}
}
