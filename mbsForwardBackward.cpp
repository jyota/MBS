// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
double mbsMvar(arma::mat x, arma::vec classes)
{
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
	
	return trace(SSCP_regression * inv(SSCP_residual));	
}

// [[Rcpp::export]]
NumericVector forwardSelection(NumericMatrix x, NumericVector classes, NumericVector selectedCols, NumericVector selectedRows, Function mvar)
{
	double maxGain = 0.0;
	//double currentT2 = mbsMvar(as<arma::mat>(x), as<arma::vec>(classes));
	int poolSize = x.ncol();
	IntegerVector pool(poolSize - selectedCols.size());
	return wrap(mbsMvar(as<arma::mat>(x), as<arma::vec>(classes)));
}
