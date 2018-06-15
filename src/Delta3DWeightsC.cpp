#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube Delta3DWeightsC(arma::cube vx,Rcpp::NumericVector Datasample) {
  /* Beware of changes:
   * This function was existing 2 times before refactoring by FP
   * and they differed by input type (NumericVector instead of cube)
   * 
   * now just the arma::cube variant was kept. This SHOULD work.
   * If problems arise: make another with original input type, 
   * rename the parameter just x and add the 2 lines
   * of code below as the function with rcpp export annotated.
   * This function stays for internal use in trainstepC
   */
  //Rcpp::IntegerVector vx_dims = x.attr("dim");
  //arma::cube vx(x.begin(), vx_dims[0], vx_dims[1], vx_dims[2], false);
  
  for (unsigned int i = 0; i < vx.n_slices; i++) {
    vx.slice(i)=vx.slice(i)-Datasample(i);
  }
  
  return vx;
}