// ------------------------------------------------------------------
// General utilities - not directly related to but used by gRips.
// ------------------------------------------------------------------

#include "RcppArmadillo.h"
#include <math.h>

using namespace Rcpp;
using namespace arma;

typedef Rcpp::NumericVector   num_vec;
typedef Rcpp::IntegerVector   int_vec;
typedef Rcpp::CharacterVector chr_vec;

//[[Rcpp::export(.c_clone)]]
SEXP clone_(SEXP& x)
{
  SEXP x2 = clone(x);
  return(x2);
}

//[[Rcpp::export]]
chr_vec list_names_(List lst){
  chr_vec out(0);
  if (lst.length() > 0)
    out = as<chr_vec>(lst.names());
  return out;
}

//[[Rcpp::export]]
mat as_emat2amat_(umat emat, int d){
  mat amat = zeros(d, d);
  uvec eids = sub2ind(size(amat), emat);
  vec vals = ones(emat.n_cols, 1);
  amat(eids) = vals;  
  amat = amat + amat.t();
  return amat;
}

//[[Rcpp::export]]
umat as_emat_complement_(umat emat, int d){

  mat amat = as_emat2amat_(emat, d);
  amat = amat - 1;
  amat = trimatl(amat);
  amat.diag().zeros();
  uvec indices = find(amat < 0);
  umat ematc   = ind2sub(size(amat), indices);
  return ematc;  
}























// ---------------------------------------------------------------------
// *** Find unique rows and columns in arma matrix ***
// ---------------------------------------------------------------------
//
// https://stackoverflow.com/questions/37143283/finding-unique-rows-in-armamat

// template <typename T>
// inline bool approx_equal_cpp(const T& lhs, const T& rhs, double tol = 0.00000001) {
//return arma::approx_equal(lhs, rhs, "absdiff", tol);
// }

// // [[Rcpp::export]]
// arma::mat unique_rows(const arma::mat& m) {

//   arma::uvec ulmt = arma::zeros<arma::uvec>(m.n_rows);

//   for (arma::uword i = 0; i < m.n_rows; i++) {
//     for (arma::uword j = i + 1; j < m.n_rows; j++) {
//       if (approx_equal_cpp(m.row(i), m.row(j))) { ulmt(j) = 1; break; }
//     }
//   }
//   return m.rows(find(ulmt == 0));
// }

// // [[Rcpp::export]]
// arma::mat unique_cols(const arma::mat& m) {

//   arma::uvec vlmt = arma::zeros<arma::uvec>(m.n_cols);
//   // Rf_PrintValue(wrap(vlmt));
//   for (arma::uword i = 0; i < m.n_cols; i++) {
//     for (arma::uword j = i + 1; j < m.n_cols; j++) {
//       if (approx_equal_cpp(m.col(i), m.col(j))) {
// 	vlmt(j) = 1; break;
//       }
//     }
//   }
//   // Rf_PrintValue(wrap(vlmt));		
//   return m.cols(find(vlmt == 0));
// }

// ---------------------------------------------------------------------


// // Returns index of first occurence of 'st' in 'x' and -1 if 'st' can
// // not be found.
// //[[Rcpp::export]]
// int find_str_ (const char* st, chr_vec x){

//   // Rcout << st << endl;
//   chr_vec::iterator it = std::find(x.begin(), x.end(), st);
//   auto out =  (it-x.begin());
//   // Rcout << out  << endl;
//   //* (i-x.begin());
//   if (it == x.end()) out=-1;
  
//   // Rcout << out  << endl;
//   return out;  
// }





// ---------------------------------------------------------------------
// *** Convert list of vectors into edge matrix (a 2xp matrix) ***
// ---------------------------------------------------------------------

// THIS FN IS VERY SLOW; Probably way faster to create amat first.

// //[[Rcpp::export]] 
// arma::mat list_to_emat (const List& lst, int shift=1){

//   int nc=lst.length(), nr=0,  m;
//   if (nc==0) return(mat(2,0));

//   nr = as<vec>(lst[0]).size();  
//   for (int i=0; i<lst.length(); ++i){
//     m = as<vec>(lst[i]).size();
//     if (m != nr)
//       stop("element %i is of different length than the previous\n", i+1);
//   }

//   mat result_mat(nr, nc);
//   size_t kk = 0;
   
//   for (int i = 0; i < lst.length(); ++i){
//     vec gg = sort(as<vec>(lst[i])) - shift;
//     result_mat.col(kk++) = gg;
//   }

//   result_mat = unique_cols(result_mat);
//   return result_mat;
// }

