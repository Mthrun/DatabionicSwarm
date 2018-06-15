using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
NumericVector lmC(NumericVector x,NumericVector yr) {
 // doku see lm() in R 
  NumericMatrix Xr(yr.length(),2);
  for(int i=0;i<yr.length();i++){
    Xr(i,0)=1;
    Xr(i,1)=x(i);
  }
  int n = Xr.nrow(), k = Xr.ncol();
  
  arma::mat X(Xr.begin(), n, k, false);       // reuses memory and avoids extra copy
  arma::colvec y(yr.begin(), yr.size(), false);
  
  arma::colvec coef = arma::solve(X, y);      // fit model y ~ X
  
  NumericVector coefs(coef.begin(),coef.end());
  
  return(coefs);
}
