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
arma::mat list2Emat_ (const List& E, int shift=1);

SEXP clone_(SEXP& x);
chr_vec list_names_(List lst);
int find_str_ (const char* st, chr_vec x);
