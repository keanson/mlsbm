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

// [[Rcpp::export]]
LogicalVector is_any(NumericVector x, double c)
{
  return x > c;
}

// [[Rcpp::export]]
double whichmax(NumericVector x)
{
  return which_max(x);
}

// [[Rcpp::export]]
NumericVector fixpi(NumericVector pi_star)
{
  NumericVector pi_star2 (pi_star.length());
  if(any(abs(pi_star) > 700).is_true())
  {
    pi_star2[which_max(pi_star)] = 1;
  }
  else
  {
    pi_star2 = exp(pi_star) / sum(exp(pi_star));
  }
  return pi_star2;
}