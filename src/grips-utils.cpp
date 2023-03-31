#include "RcppArmadillo.h"
#include "general-utils.h"

using namespace Rcpp;
using namespace arma;

typedef Rcpp::NumericVector   num_vec;
typedef Rcpp::IntegerVector   int_vec;
typedef Rcpp::CharacterVector chr_vec;

// ---------------------------------------------------------------------
// *** gRips utilities ***
// ---------------------------------------------------------------------


// [[Rcpp::export]]
mat project_K_onto_G_(mat& K, umat& emc){

  uvec r0 = {0}, r1 = {1};
  mat emc2 = conv_to<mat>::from(emc);
  uvec s0 = conv_to<uvec>::from(emc2.rows(r0));
  uvec s1 = conv_to<uvec>::from(emc2.rows(r1));
  for (size_t j=0; j<emc.n_cols; j++){
    // Rcout << s0[j] << "  " << s1[j] << "\n"; 
    K(s0[j], s1[j]) = 0;
    K(s1[j], s0[j]) = 0;
  }
  return(K);
}

// [[Rcpp::export]]
double mnormone_(mat& Delta){
  rowvec s = sum(abs(Delta));
  // s.print();
  return(max(s));
}
// mnormone(Delta) = max_v(\sum_u |Delta_{uv}|)


double get_conv_ref(const List& aux){
  CharacterVector vn = list_names_(aux);
  if (vn.length()){
    if (find_str_("conv_ref", vn) < 0)
      stop("'conv_ref' not found in 'aux'");
  }
  double out = as<double>(aux["conv_ref"]);
  return out;
}



mat initSigma_(mat& S)
{
  int p = S.n_rows;
  vec s = S.diag();
  mat Sigma(p, p, fill::eye);
  Sigma.diag() = s;
  return Sigma;
}

mat initK_(mat& S)
{
  int p = S.n_rows;
  vec s = S.diag();
  mat K(p, p, fill::eye);
  K.diag() = 1/s;
  return K;
}


int method2int_(CharacterVector method){
  CharacterVector options =
    CharacterVector::create("ips", "mtp2", "lasso", "hybrid");
  
  IntegerVector m =  match(method, options) - 1;
  int v = (int) m[0];
  return v;
}




// [[Rcpp::export]]
double ggm_logL_(mat& S, mat& K, int nobs)
{
  double trKS = accu(K % S);
  // Rf_PrintValue(wrap(trKS));
  double val, sign;
  log_det(val, sign, K);
  // double logL = nobs * (log(det(K)) - trKS) / 2;
  double logL = nobs * (val - trKS) / 2;

  return logL;
}


