// ------------------------------------------------------------------
// General utilities - not directly related to but used by gRips.
// ------------------------------------------------------------------

#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;

typedef Rcpp::NumericVector   num_vec;
typedef Rcpp::IntegerVector   int_vec;
typedef Rcpp::CharacterVector chr_vec;

arma::mat unique_rows(const arma::mat& m);
arma::mat unique_cols(const arma::mat& m);
arma::mat list_to_emat (const List& lst, int shift=1);
mat as_emat2amat_(umat emat, int d);
umat as_emat_complement_(umat emat, int d);


SEXP clone_(SEXP& x);
chr_vec list_names_(List lst);
int find_str_ (const char* st, chr_vec x);
