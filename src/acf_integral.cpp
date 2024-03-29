#include <Rcpp.h>
using namespace Rcpp;

//' Multiplies by two!
//'
//' @param x number
//' @return twice as much as `x`
//'
//' @export
// [[Rcpp::export]]
int times_two(int x) {
  return x * 2.0;
}

double main_cross_moment(float x, float y, float beta, float b_oh, float b_o_exc_h) {
  double result(pow(1. + (x+y)/beta, b_oh));
  result *= pow(1. + x/beta, b_o_exc_h);
  result *= pow(1. + y/beta, b_o_exc_h);

  return result;
}

double main_moment(float x, float beta, float alpha) {
  return pow(1. + x / beta, -alpha);
}

//' Computes cross moment in trawl gamma-exponential mixture
//'
//' @param xs Grid to integrate on.
//' @param delta Grid mesh size
//' @param beta GPD scale parameter in the sense of Noven et al.
//' @param b_oh - alpha * mu(A inter A_h) / mu  (A)
//' @param b_o_exc_h - alpha * mu(A excl. A_h) / mu  (A)
//' @return Cross moment value
//'
//' @export
// [[Rcpp::export]]
double cross_moment(NumericVector xs, float delta, float beta, float b_oh, float b_o_exc_h) {
  double total = 0;
  NumericVector ys = xs + delta / 2.;

  for(NumericVector::iterator i = ys.begin(); i != ys.end(); ++i) {
    for(NumericVector::iterator j = ys.begin(); j != ys.end(); ++j) {
      total += main_cross_moment(*i, *j, beta, b_oh, b_o_exc_h) * pow(delta, 2.);
    }
  }
  return total;
}

//' Computes first moment in trawl gamma-exponential mixture
//'
//' @param xs Grid to integrate on.
//' @param delta Grid mesh size
//' @param beta GPD scale parameter in the sense of Noven et al.
//' @param b_oh - alpha * mu(A inter A_h) / mu  (A)
//' @param b_o_exc_h - alpha * mu(A excl. A_h) / mu  (A)
//' @return First moment value
//'
//' @export
// [[Rcpp::export]]
double first_moment(NumericVector xs, float delta, float beta, float b_oh, float b_o_exc_h) {
  double total = 0, alpha(-b_oh-b_o_exc_h);
  NumericVector ys = xs + delta / 2.;
  for(NumericVector::iterator i = ys.begin(); i != ys.end(); ++i) {
    total += main_moment(*i, beta, alpha) * delta;
  }
  return total;
}

//' Computes square moment in trawl gamma-exponential mixture
//'
//' @param xs Grid to integrate on.
//' @param delta Grid mesh size
//' @param beta GPD scale parameter in the sense of Noven et al.
//' @param b_oh - alpha * mu(A inter A_h) / mu  (A)
//' @param b_o_exc_h - alpha * mu(A excl. A_h) / mu  (A)
//' @return Square moment value
//'
//' @export
// [[Rcpp::export]]
double square_moment(NumericVector xs, float delta, float beta, float b_oh, float b_o_exc_h) {
  double total = 0, alpha(-b_oh-b_o_exc_h);
  NumericVector ys = xs + delta / 2.;
  for(NumericVector::iterator i = ys.begin(); i != ys.end(); ++i) {
    for(NumericVector::iterator j = ys.begin(); j != ys.end(); ++j) {
      total += main_cross_moment(*i, *j, beta, -alpha, 0.0) * pow(delta, 2.);
    }
  }
  return total;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
microbenchmark::microbenchmark(times_two(42))
microbenchmark::microbenchmark(cross_moment(c(1,2,3,4), 2., .1, 1., -.3))
microbenchmark::microbenchmark(first_moment(c(1,2,3,4), 2., .1, 1., -.3))
*/
