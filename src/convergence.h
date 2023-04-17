#include "RcppArmadillo.h"
using namespace Rcpp;

double max_abs_(const mat& S);
double max_abs_diag_(const mat& S);
double max_abs_diff_(const mat& S, const mat& Sigma);
double max_abs_diff_rel_(const mat& S, const mat& Sigma);
double max_abs_diag_diff_(const mat& S, const mat& Sigma);
double max_diag_diff_(const mat& S, const mat& Sigma);

double max_abs_diff_on_emat_  (const mat& S, const mat& Sigma, const umat& E, int shift=1);
double max_abs_diff_on_elst_ (const mat& S, const mat& Sigma, const List& E, int shift=1);

double mean_abs_diff_on_emat_ (const mat& S, const mat& Sigma, const umat& E, int shift=1);
double mean_abs_diff_on_elst_(const mat& S, const mat& Sigma, const List& E, int shift=1);

