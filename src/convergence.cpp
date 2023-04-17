#include "RcppArmadillo.h"
#include "grips_utils.h"
#include "general_utils.h"
#include <math.h>

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::interfaces(r,cpp)]]

using namespace Rcpp;
using namespace arma;

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)>(b))?(b):(a))



// --- Utilities ---
// parallel functions exist in R

//[[Rcpp::export]] 
double max_abs_(const mat& S)
{
  double d = max(max(abs(S)));
  return d;
}

//[[Rcpp::export]] 
double max_abs_diag_(const mat& S)
{
  double d = as_scalar(max(abs(S.diag())));
  return d;
}

//[[Rcpp::export]] 
double max_abs_diff_rel_(const mat& S, const mat& Sigma)
{
  double d = max(max(abs(S - Sigma)));
  double r = max(max(abs(Sigma)));
  return d / r;
}

//[[Rcpp::export]] 
double max_abs_diff_(const mat& S, const mat& Sigma)
{
  // double d = max(max(abs(S - Sigma)));
  // return d;

  double d=0;
  for (size_t i=0; i<S.n_rows; i++){
    for (size_t j=0; j<=i; j++){
      d = MAX(d, fabs(S(i,j) - Sigma(i,j)));    
    }
  }
  return d;
}

//[[Rcpp::export]] 
double max_abs_diag_diff_(const mat& S, const mat& Sigma)
{
  double d=0;
  for (size_t i=0; i<S.n_rows; i++){
    d = MAX(d, fabs(S(i,i) - Sigma(i,i)));    
  }
  return d;  
 
  // vec d = S.diag();
  // d -= Sigma.diag();
  // return max(abs(d));
}

//[[Rcpp::export]] 
double max_diag_diff_(const mat& S, const mat& Sigma)
{
  double d=0;
  for (size_t i=0; i<S.n_rows; i++){
    d = MAX(d, S(i,i) - Sigma(i,i));    
  }
  return d;  
  
  // vec d = S.diag();
  // d -= Sigma.diag();
  // return max(d);
}



// --- END Utilities ---


// ----------------------------------------------------------------------
// Implements ||S(G) - Sigma(G)||_\infty
//
// Likelihood equation condition for GGM
// ----------------------------------------------------------------------


//[[Rcpp::export]] 
vec diff_on_emat_(const mat& S, const mat& Sigma, const umat& E, int shift=1){
  uvec eids = sub2ind(size(S), E - shift);           // Obtain Element IDs
  vec v = S.elem(eids) - Sigma.elem(eids);   // Values of the Elements
  vec d = S.diag() - Sigma.diag();
  vec out = join_cols(d, v);
  return out;
}

//[[Rcpp::export]] 
vec diff_on_elst_(const mat& S, const mat& Sigma, const List& E, int shift=1){
  umat emat = conv_to<umat>::from(list_to_emat(E, shift));
  // emat.print();
  return diff_on_emat_(S, Sigma, emat, 0);
  // return diff_on_emat_(S, Sigma, emat, shift);

  // uvec eids = sub2ind( size(S), emat);           // Obtain Element IDs
  // vec v = S.elem( eids ) - Sigma.elem( eids );   // Values of the Elements
  // vec d = S.diag() - Sigma.diag();
  // vec out = join_cols(d, v);
  // return out;
  
}

//[[Rcpp::export]] 
double max_abs_diff_on_emat_(const mat& Sigma, const mat& S, const umat& E, int shift=1){
  vec d = diff_on_emat_(Sigma, S, E, shift);
  return max(abs(d));
}

//[[Rcpp::export]] 
double mean_abs_diff_on_emat_(const mat& Sigma, const mat& S, const umat& E, int shift=1){
  vec d = diff_on_emat_(Sigma, S, E, shift);
  return (sum(abs(d)) / d.size());
}

//[[Rcpp::export]] 
double max_abs_diff_on_elst_(const mat& Sigma, const mat& S, const List& E, int shift=1){
  vec d = diff_on_elst_(Sigma, S, E, shift);
  return max(abs(d));
}

//[[Rcpp::export]] 
double mean_abs_diff_on_elst_(const mat& Sigma, const mat& S, const List& E, int shift=1){
  vec d = diff_on_elst_(Sigma, S, E, shift);
  return (sum(abs(d)) / d.size());
}

//[[Rcpp::export]] 
double max_diff_on_emat_(const mat& Sigma, const mat& S, const mat& E){
  double out=0, dd; //, v1, v2;
  int u, v;

  out = max_diag_diff_(Sigma, S);

  for (size_t j=0; j<E.n_cols; ++j){
    u = E(0, j) - 1; v = E(1, j) - 1 ;
    // v1 = Sigma(u, v);   v2 = S(u, v);    dd  = (v1 - v2);
    dd  = Sigma(u, v) - S(u, v);
    // Rcout << v1 << " " << v2 << " " << out << std::endl;
    out = MAX(out, dd);	    
  }
  return out;
}	


// ----------------------------------------------------------------------
// implements ||{S(G) - Sigma(G)}^+||_\infty
//
// Likelihood equation condition for MTP2
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
// Implements || { S(G) - Sigma(G) } * K(G) ||_\infty
//
// Likelihood equation condition for MTP2
// ----------------------------------------------------------------------

// //[[Rcpp::export]] 
// double max_abs_diff_on_EK_(const mat& S, const mat& Sigma, const mat& E, const mat& K){
//   double out=0, dd; 
//   int u, v;
  
//   for (size_t i=0; i<Sigma.n_rows; ++i)
//     out = MAX(out, fabs(Sigma(i, i) - S(i, i)) * K(i, i));
 
//   for (size_t j=0; j<E.n_cols; ++j){
//     u = E(0,j) - 1; v = E(1,j) - 1 ;
//     // v1 = Sigma(u, v);    v2 = S(u, v);    v3 = K(u, v);    dd  = fabs((v1 - v2) * v3);
//     dd  = fabs((Sigma(u, v) - S(u, v)) * K(u, v));
//     //Rcout << v1 << " " << v2 << " " << out << std::endl;
//     out = MAX(out, dd);	    
//   }
//   return out;
// }	






// //[[Rcpp::export]] 
// double califa_(const mat& S, const mat& Sigma, const mat& Sigmaold)
// {
//   double N = S.n_rows;
//   double d = 0; //accu(abs(Sigma - Sigmaold)); // / (N * N);
//   //double r = max(max(abs(Sigma)));
//   double r = 0;
//   for (int i=0; i<N; i++){
//     for (int j=0; j<=i; j++){
//       // double w = fabs(S(i,j));
//       if (i != j)
// 	r += fabs(S(i,j));

//       d += fabs(Sigma(i,j) - Sigmaold(i,j));
//       // Rcout <<  " i : " << i << " j : " << j <<  " w ::: " << w << "  r ::: " << r << std::endl;
//     }
//   }
//   r = r / (N * (N-1) / 2);
//   d = d / (N * (N+1) / 2);
//   double rat = d / r;

//   // Rcout << "d : " << d << " r : " << r << " rat : " << rat << std::endl;
  
//   return rat;
// }
