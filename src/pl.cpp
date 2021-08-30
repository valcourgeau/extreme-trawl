#include <cstdio>
#include <cmath>
#include <Rcpp.h>

using namespace Rcpp;

const double EPSILON = std::numeric_limits<double>::epsilon();

//' Cpp pairwise likelihood 0-0 (zero-zero).
//'
//' @param alpha Gamma shape.
//' @param beta Gamma rate.
//' @param kappa Exceedance controlling parameter.
//' @param b_1 B one.
//' @param b_2 B two.
//' @param b_3 B three.
//' @return pairwise likelihood 0-0 (zero-zero).
// [[Rcpp::export]]
double cpp_case_zero_zero(double alpha, double beta, double kappa, double b_1, double b_2, double b_3) {
  double A(b_1+b_3);
  double b1(-alpha*b_1/A), b2(-alpha*b_2/A), b3(-alpha*b_3/A);
  double tmp(1.0 - 2.0*pow(1+kappa/beta, -alpha) + pow(1.0+kappa/beta, b1+b3)*pow(1.0+2.0*kappa/beta, b2));

  return tmp;
}

//' Cpp pairwise likelihood 1-0 (one-zero).
//'
//' @param xs Data.
//' @param alpha Gamma shape.
//' @param beta Gamma rate.
//' @param kappa Exceedance controlling parameter.
//' @param b_1 B one.
//' @param b_2 B two.
//' @param b_3 B three.
//' @return pairwise likelihood 1-0 (one-zero).
// [[Rcpp::export]]
double cpp_case_one_zero(NumericVector xs, double alpha, double beta, double kappa, double b_1, double b_2, double b_3) {
  double A(b_1+b_3);
  double b1(-alpha*b_1/A), b2(-alpha*b_2/A);
  double tmp(alpha/beta * pow(1.0+kappa/beta, -alpha-1.0));

  if (std::abs(xs[0]) < EPSILON){
    xs[0] = xs[1];
    xs[1] = 0.0;
  }
  double x(xs[0]);

  tmp += 1.0/beta*pow(1.0+(kappa+x)/beta, b1-1.0)*pow(1.0+(2.0*kappa+x)/beta, b2-1.0)*pow(1.0+kappa/beta, b1)*(-alpha*(1.0+(kappa+x)/beta)+b1*kappa/beta);

  return tmp;
}

//' Cpp pairwise likelihood 1-1 (one-one).
//'
//' @param xs Data.
//' @param alpha Gamma shape.
//' @param beta Gamma rate.
//' @param kappa Exceedance controlling parameter.
//' @param b_1 B one.
//' @param b_2 B two.
//' @param b_3 B three.
//' @return pairwise likelihood 1-1 (one-one).
// [[Rcpp::export]]
double cpp_case_one_one(NumericVector xs, double alpha, double beta, double kappa, double b_1, double b_2, double b_3) {
  double A(b_1+b_3);
  double b1(-alpha*b_1/A), b2(-alpha*b_2/A), b3(-alpha*b_3/A);
  double x1(xs[1]), x2(xs[2]);
  double tmp_1(1.0/pow(beta, 2.0)*pow(1.0+(2.0*kappa+x1+x2)/beta, b2-2.0)*pow(1.0+(kappa+x1)/beta, b1-1.0)*pow(1+(kappa+x2)/beta, b3-1.0));
  double tmp_2(b2*(b2-1.0)*(1.0+(kappa+x1)/beta)*(1.0+(kappa+x2)/beta)+pow(b1, 2.0)*pow(1.0+(2.0*kappa+x1+x2)/beta, 2.0));
  double tmp_3(b2*b1*(1.0+(2*kappa+x1+x2)/beta)*((1.0+(kappa+x1)/beta)+(1.0+(kappa+x2)/beta)));

  return tmp_1 * (tmp_2 + tmp_3);
}

//' Cpp pairwise likelihood separator that chooses which routine to use.
//'
//' @param xs Data.
//' @param alpha Gamma shape.
//' @param beta Gamma rate.
//' @param kappa Exceedance controlling parameter.
//' @param b_1 B one.
//' @param b_2 B two.
//' @param b_3 B three.
//' @return PL likelihood values for `xs`.
// [[Rcpp::export]]
double cpp_case_separator(NumericVector xs, double alpha, double beta, double kappa, double b_1, double b_2, double b_3) {
  assert(xs.size() == 2);
  double product = std::accumulate(xs.begin(), xs.end(), 1, std::multiplies<double>());
  bool all_nonpositive = std::all_of(xs.begin(), xs.end(), [](double x){ return x <= EPSILON; });

  if(all_nonpositive){
    return cpp_case_zero_zero(alpha, beta, kappa, b_1, b_2, b_3);
  }else if(std::abs(product) < EPSILON){
    return cpp_case_one_zero(xs, alpha, beta, kappa, b_1, b_2, b_3);
  }else{
    return cpp_case_one_one(xs, alpha, beta, kappa, b_1, b_2, b_3);
  }
}
