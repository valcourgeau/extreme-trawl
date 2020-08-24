#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int timesTwo(int x) {
  return x * 2;
}

double main_cross_moment(float x, float y,float beta, float b_oh, float b_o_exc_h) {
  double result(pow(1. + (x+y)/beta, b_oh));
  result *= pow(1. + x/beta, b_o_exc_h);
  result *= pow(1. + y/beta, b_o_exc_h);

  return result;
}

double main_moment(float x, float beta, float alpha) {

  return pow(1. + x / beta, -alpha);
}

// [[Rcpp::export]]
double CrossMoment(NumericVector xs, float delta, float beta, float b_oh, float b_o_exc_h) {
  double total = 0;
  NumericVector ys = xs + delta / 2.;

  for(NumericVector::iterator i = ys.begin(); i != ys.end(); ++i) {
    for(NumericVector::iterator j = ys.begin(); j != ys.end(); ++j) {
      total += main_cross_moment(*i, *j, beta, b_oh, b_o_exc_h) * pow(delta, 2.);
    }
  }
  return total;
}

// [[Rcpp::export]]
double FirstMoment(NumericVector xs, float delta, float beta, float b_oh, float b_o_exc_h) {
  double total = 0, alpha(-b_oh-b_o_exc_h);
  for(NumericVector::iterator i = xs.begin(); i != xs.end(); ++i) {
    total += main_moment(*i, beta, alpha) * delta;
  }
  return total;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
main_cross_moment(1., 2., .1, 1., -.3, -.4)
*/
