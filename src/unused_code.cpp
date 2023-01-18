#include "RcppArmadillo.h"
#include "grips-utils.h"
#include <math.h>

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::interfaces(r,cpp)]]

using namespace Rcpp;
using namespace arma;


// ---------------------------------------------------------------------
// *** Converts list of vectors to a matrix ***
// ---------------------------------------------------------------------
//
// Either column by column or row by row.  All vectors must have same
// length. Taken from SO, but I forgot to copy url.

// [[Rcpp::export]]
Rcpp::NumericMatrix list2row_(Rcpp::List input_list, int d=2){

  unsigned int n = input_list.length();
  
  if (n == 0){
    Rcpp::NumericMatrix result_mat = Rcpp::no_init(0, d);    
    return  result_mat;
  }
  
  Rcpp::NumericVector testvals = input_list[0];
  unsigned int elems = testvals.length();
  
  Rcpp::NumericMatrix result_mat = Rcpp::no_init(n, elems);
  
  // fill by row
  for(unsigned int i = 0; i < n; i++) {
    Rcpp::NumericVector row_val = input_list[i];
    
    if(elems != row_val.length()) {
      Rcpp::stop("Length of row does not match matrix requirements"); 
    }
    result_mat(i, Rcpp::_) = row_val;
  }
  
  return result_mat;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix list2col_(Rcpp::List input_list, int d=2){
  return transpose(list2row_(input_list, d));
}

// ---------------------------------------------------------------------

