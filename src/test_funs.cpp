// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double POW(NumericMatrix Ps, double p) {
  return pow(Ps(0,0),p);
}

// [[Rcpp::export]]
double logPow(NumericMatrix Ps, double p) {
  return log(pow(Ps(0,0),p));
}

// [[Rcpp::export]]
LogicalVector isFinite(NumericVector x)
{
  // return arma::is_finite(x);
  return is_na(x);
}

// [[Rcpp::export]]
NumericVector fix_NAs(NumericVector x)
{
  x[is_na(x)] = 1;
  return x;
}