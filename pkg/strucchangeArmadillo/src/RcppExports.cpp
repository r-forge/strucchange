#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// sc_cpp_recresid
NumericVector sc_cpp_recresid(const arma::mat& X, const arma::vec& y, unsigned int start, unsigned int end, const double& tol, const double& rcond_min);
RcppExport SEXP strucchange_sc_cpp_recresid(SEXP XSEXP, SEXP ySEXP, SEXP startSEXP, SEXP endSEXP, SEXP tolSEXP, SEXP rcond_minSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< unsigned int >::type start(startSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type end(endSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const double& >::type rcond_min(rcond_minSEXP);
    rcpp_result_gen = Rcpp::wrap(sc_cpp_recresid(X, y, start, end, tol, rcond_min));
    return rcpp_result_gen;
END_RCPP
}

