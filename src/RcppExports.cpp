// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// times_two
int times_two(int x);
RcppExport SEXP _extreme_trawl_times_two(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(times_two(x));
    return rcpp_result_gen;
END_RCPP
}
// cross_moment
double cross_moment(NumericVector xs, float delta, float beta, float b_oh, float b_o_exc_h);
RcppExport SEXP _extreme_trawl_cross_moment(SEXP xsSEXP, SEXP deltaSEXP, SEXP betaSEXP, SEXP b_ohSEXP, SEXP b_o_exc_hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xs(xsSEXP);
    Rcpp::traits::input_parameter< float >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< float >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< float >::type b_oh(b_ohSEXP);
    Rcpp::traits::input_parameter< float >::type b_o_exc_h(b_o_exc_hSEXP);
    rcpp_result_gen = Rcpp::wrap(cross_moment(xs, delta, beta, b_oh, b_o_exc_h));
    return rcpp_result_gen;
END_RCPP
}
// first_moment
double first_moment(NumericVector xs, float delta, float beta, float b_oh, float b_o_exc_h);
RcppExport SEXP _extreme_trawl_first_moment(SEXP xsSEXP, SEXP deltaSEXP, SEXP betaSEXP, SEXP b_ohSEXP, SEXP b_o_exc_hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xs(xsSEXP);
    Rcpp::traits::input_parameter< float >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< float >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< float >::type b_oh(b_ohSEXP);
    Rcpp::traits::input_parameter< float >::type b_o_exc_h(b_o_exc_hSEXP);
    rcpp_result_gen = Rcpp::wrap(first_moment(xs, delta, beta, b_oh, b_o_exc_h));
    return rcpp_result_gen;
END_RCPP
}
// square_moment
double square_moment(NumericVector xs, float delta, float beta, float b_oh, float b_o_exc_h);
RcppExport SEXP _extreme_trawl_square_moment(SEXP xsSEXP, SEXP deltaSEXP, SEXP betaSEXP, SEXP b_ohSEXP, SEXP b_o_exc_hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xs(xsSEXP);
    Rcpp::traits::input_parameter< float >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< float >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< float >::type b_oh(b_ohSEXP);
    Rcpp::traits::input_parameter< float >::type b_o_exc_h(b_o_exc_hSEXP);
    rcpp_result_gen = Rcpp::wrap(square_moment(xs, delta, beta, b_oh, b_o_exc_h));
    return rcpp_result_gen;
END_RCPP
}
// cpp_case_zero_zero
double cpp_case_zero_zero(double alpha, double beta, double kappa, double b_1, double b_2, double b_3);
RcppExport SEXP _extreme_trawl_cpp_case_zero_zero(SEXP alphaSEXP, SEXP betaSEXP, SEXP kappaSEXP, SEXP b_1SEXP, SEXP b_2SEXP, SEXP b_3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type b_1(b_1SEXP);
    Rcpp::traits::input_parameter< double >::type b_2(b_2SEXP);
    Rcpp::traits::input_parameter< double >::type b_3(b_3SEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_case_zero_zero(alpha, beta, kappa, b_1, b_2, b_3));
    return rcpp_result_gen;
END_RCPP
}
// cpp_case_one_zero
double cpp_case_one_zero(NumericVector xs, double alpha, double beta, double kappa, double b_1, double b_2, double b_3);
RcppExport SEXP _extreme_trawl_cpp_case_one_zero(SEXP xsSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP kappaSEXP, SEXP b_1SEXP, SEXP b_2SEXP, SEXP b_3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xs(xsSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type b_1(b_1SEXP);
    Rcpp::traits::input_parameter< double >::type b_2(b_2SEXP);
    Rcpp::traits::input_parameter< double >::type b_3(b_3SEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_case_one_zero(xs, alpha, beta, kappa, b_1, b_2, b_3));
    return rcpp_result_gen;
END_RCPP
}
// cpp_case_one_one
double cpp_case_one_one(NumericVector xs, double alpha, double beta, double kappa, double b_1, double b_2, double b_3);
RcppExport SEXP _extreme_trawl_cpp_case_one_one(SEXP xsSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP kappaSEXP, SEXP b_1SEXP, SEXP b_2SEXP, SEXP b_3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xs(xsSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type b_1(b_1SEXP);
    Rcpp::traits::input_parameter< double >::type b_2(b_2SEXP);
    Rcpp::traits::input_parameter< double >::type b_3(b_3SEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_case_one_one(xs, alpha, beta, kappa, b_1, b_2, b_3));
    return rcpp_result_gen;
END_RCPP
}
// cpp_case_separator
double cpp_case_separator(NumericVector xs, double alpha, double beta, double kappa, double b_1, double b_2, double b_3);
RcppExport SEXP _extreme_trawl_cpp_case_separator(SEXP xsSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP kappaSEXP, SEXP b_1SEXP, SEXP b_2SEXP, SEXP b_3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xs(xsSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< double >::type b_1(b_1SEXP);
    Rcpp::traits::input_parameter< double >::type b_2(b_2SEXP);
    Rcpp::traits::input_parameter< double >::type b_3(b_3SEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_case_separator(xs, alpha, beta, kappa, b_1, b_2, b_3));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_extreme_trawl_times_two", (DL_FUNC) &_extreme_trawl_times_two, 1},
    {"_extreme_trawl_cross_moment", (DL_FUNC) &_extreme_trawl_cross_moment, 5},
    {"_extreme_trawl_first_moment", (DL_FUNC) &_extreme_trawl_first_moment, 5},
    {"_extreme_trawl_square_moment", (DL_FUNC) &_extreme_trawl_square_moment, 5},
    {"_extreme_trawl_cpp_case_zero_zero", (DL_FUNC) &_extreme_trawl_cpp_case_zero_zero, 6},
    {"_extreme_trawl_cpp_case_one_zero", (DL_FUNC) &_extreme_trawl_cpp_case_one_zero, 7},
    {"_extreme_trawl_cpp_case_one_one", (DL_FUNC) &_extreme_trawl_cpp_case_one_one, 7},
    {"_extreme_trawl_cpp_case_separator", (DL_FUNC) &_extreme_trawl_cpp_case_separator, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_extreme_trawl(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
