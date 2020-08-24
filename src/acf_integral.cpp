#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x, float delta) {
  return x * 2;
}

double main_cross_moment(float x, float y, float delta, float beta, float b_oh, float b_o_exc_h) {
  x += delta / 2.;
  y += delta / 2.;
  double result(pow(1. + (x+y)/beta, b_oh+b_o_exc_h));
  result *= pow(1. + x/beta, b_o_exc_h);
  result *= pow(1. + y/beta, b_o_exc_h);
  return result;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
main_cross_moment(1., 2., .1, 1., -.3, -.4)
*/
