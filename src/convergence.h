#include "RcppArmadillo.h"
using namespace Rcpp;

double max_abs_(const mat& S);
double max_abs_diag_(const mat& S);
double max_abs_diff_(const mat& S, const mat& Sigma);
double max_abs_diff_rel_(const mat& S, const mat& Sigma);
double max_abs_diag_diff_(const mat& S, const mat& Sigma);
double max_diag_diff_(const mat& S, const mat& Sigma);

bool does_edge_fit_to_data (mat& S, uvec& uv, mat& K, mat& Sigma, double& eps, mat& lambda, int imeth=0);
bool does_model_fit_to_data(mat& S, umat& EE,  mat& K, mat& Sigma, double& eps, mat& lambda, int imeth=0);
bool does_model_fit_to_data_relaxed(mat& S, mat& Sigma, mat& Sigma_prev, double& eps);

double max_abs_diff_on_emat_  (const mat& S, const mat& Sigma, const umat& E, int shift=1);
double max_abs_diff_on_Elist_ (const mat& S, const mat& Sigma, const List& E, int shift=1);

double mean_abs_diff_on_emat_ (const mat& S, const mat& Sigma, const umat& E, int shift=1);
double mean_abs_diff_on_Elist_(const mat& S, const mat& Sigma, const List& E, int shift=1);

double califa_(const mat& S, const mat& Sigma, const mat& Sigmaold);



