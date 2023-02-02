#include "RcppArmadillo.h"
#include "grips-utils.h"
#include "convergence.h"
#include "arma_utils.h"
#include "precision.h"
#include <iostream>
#include <vector>
#include <algorithm> // for std::nth_element
#include <stdio.h>

using namespace Rcpp;
using namespace arma;

typedef Rcpp::NumericVector   num_vec;
typedef Rcpp::IntegerVector   int_vec;
typedef Rcpp::CharacterVector chr_vec;

//[[Rcpp::export]]
mat as_emat2amat_(umat emat, int d){
  mat amat = zeros(d, d);
  arma::uvec eids = sub2ind(size(amat), emat);
  vec vals = ones(emat.n_cols, 1);
  amat(eids) = vals;
  
  amat = amat + amat.t();
  return amat;
}

//[[Rcpp::export]]
umat as_emat_complement_(umat emat, int d, int shift=1){

  mat amat = as_emat2amat_(emat, d);
  amat = amat - 1;
  amat = trimatl(amat);
  amat.diag().zeros();
  // amat.print();

  uvec indices = find(amat < 0);
  umat ematc   = ind2sub(size(amat), indices);
  return ematc;  
}

//[[Rcpp::export]]
mat inv_qr_(mat& X){
  arma::mat Q, R;
  arma::qr(Q,R,X);
  mat Ri = inv(R);
  mat out = Ri * Q.t();
  return(out);
}

// ### CALLING R FUNCTIONS


double mean_abs_diff_non_edge_(mat& Sigma1, mat& Sigma2, umat& emc){
  double out = accu(abs(Sigma1 - Sigma2)) / (2 * emc.n_cols);
  return out;
}

// double mean_abs_diff_non_edge_R(mat& Sigma1, mat& Sigma2, umat& emc){
  // Environment gRips_env = Environment::namespace_env("gRips");
  // Function R_fun = gRips_env["mean_abs_diff_non_edge"];
  // NumericVector result = R_fun(wrap(Sigma1), wrap(Sigma2), wrap(emc+1));
  // double out=result(0);
  // return out;  
// }

double diff_fun_R(mat& Sigma, mat& K, umat emc){
  Environment gRips_env = Environment::namespace_env("gRips");
  Function R_fun = gRips_env["diff_fun"];
  NumericVector result = R_fun(wrap(Sigma), wrap(K), wrap(emc+1));
  double out=result(0);
  return out;
}


// [[Rcpp::export]]
double diff_fun_(mat& Sigma, mat& K, umat emc){
  // emc.print();
  uvec ind = sub2ind(size(Sigma), emc);
  vec Kuv = K(ind);
  // Kuv.print();
  vec diagS = Sigma.diag();

  // Rprintf("diagSigma\n");  diagS.print();

  uvec r0 = {0}, r1 = {1};
  mat emc2 = conv_to<mat>::from(emc);
  // emc2.print();
  uvec s0 = conv_to<uvec>::from(emc2.rows(r0));
  uvec s1 = conv_to<uvec>::from(emc2.rows(r1));
  // Rprintf("s0\n");  s0.print();
  // Rprintf("s1\n");  s1.print();
  vec Suu = diagS.elem(s0);
  vec Svv = diagS.elem(s1);
  // Rprintf("Suu\n");  Suu.print();
  // Rprintf("Svv\n");  Svv.print();
  int d = Sigma.n_rows * (Sigma.n_rows + 1) / 2;
  double out = sum(abs(abs(Kuv) % sqrt(Suu % Svv))) / d;
  return out;
}





// ### 

void update_row_Sigma2_(int u, mat& Sigma, const mat& amat, bool print=false){

  uvec u_    = {(unsigned int) u};      // convert int to uvec
  uvec ne_u_ = find(amat.rows(u_) > 0); // Returns column vector
  // Rcout << "ne_u_: " <<  "\n"; ne_u_.t().print();
  
  vec w(Sigma.n_cols, fill::zeros);  

  if (ne_u_.n_rows > 0){  
    mat s12 = Sigma.cols(u_);
    // vec tt = solve(Sigma.submat(ne_u_, ne_u_), s12.rows(ne_u_)); 

    mat AA = Sigma.submat(ne_u_, ne_u_);
    vec bb = s12.rows(ne_u_);
    vec tt = solve(AA, bb);
    // Rprintf("AA:\n"); AA.print();
    // Rprintf("bb:\n"); bb.t().print();
    // Rprintf("tt:\n"); tt.t().print();
    
    // vec beta(Sigma.n_cols, fill::zeros);
    // beta.elem(ne_u_) = tt;
    // w = Sigma * beta;
    // Rprintf("beta:\n"); beta.t().print();
    w = Sigma.cols(ne_u_) * tt;
  }

  // Rprintf("inserting:\n");
  double sigma_uu = Sigma(u, u); // Store this because element is overwritten below
    
  Sigma.col(u) = w;
  Sigma.row(u) = w.t();
  Sigma(u, u)  = sigma_uu;        // Restore element
  // Rprintf("inserting - done:\n");
}


void update_row_K2_(int u, mat& Sigma, mat& K, bool update_K=false, bool print=false){
  
  if (update_K){
    // int shift=0; // u is already zero-based
    // Rprintf("smart update of K\n"); K.print();
    
    ivec u__ = {(int) u};
    uvec u_  = {(unsigned int) u};
    // Rprintf("u__\n"); u__.print();

    arma::uvec uc_  = arma::linspace<arma::uvec>(0, K.n_cols - 1, K.n_cols);
    uc_.shed_rows(u_);

    // mat K_uc_u     = extract_uv_(K, -u__,  u__, shift=shift);
    // mat K_uc_uc    = extract_uv_(K, -u__, -u__, shift=shift);
    // mat Sigma_uc_u = extract_uv_(Sigma, -u__, u__, shift=shift);

    double k_uu     = as_scalar(K(u, u));    
    double sigma_uu = as_scalar(Sigma(u, u));

    mat K_uc_u  = K.submat(uc_, u_);
    mat K_uc_uc = K.submat(uc_, uc_);
    mat Sigma_uc_u = Sigma.submat(uc_, u_);
    
    // Rprintf("K_uc_u:\n"); K_uc_u.t().print();
    // Rprintf("K_uc_uc:\n"); K_uc_uc.print();
    // Rprintf("Sigma_uc_u:\n"); Sigma_uc_u.t().print();
    
    mat CC2 = K_uc_uc - K_uc_u * trans(K_uc_u) / k_uu;
    mat DD2 = CC2 * Sigma_uc_u;

    // Rprintf("CC2:\n"); CC2.print();
    // Rprintf("DD2:\n"); DD2.t().print(); 

    // DO UPDATE
    mat k_uu_upd2  = 1 / (sigma_uu - trans(Sigma_uc_u) * DD2);
    mat K_uc_u_upd = as_scalar(k_uu_upd2) * DD2;
    
    mat new2 = K_uc_u_upd * trans(K_uc_u_upd) / as_scalar(k_uu_upd2);
    mat old2 = K_uc_u     * trans(K_uc_u)     / as_scalar(k_uu);
    mat RR3  = K_uc_uc + new2 - old2;

    // arma::uvec vv_  = arma::linspace<arma::uvec>(0, K.n_cols - 1, K.n_cols);
    // arma::uvec uc_  = setdiff_(vv_, u_);
    // arma::uvec rr = vv_;
    // rr.shed_rows(u_);
    // Rprintf("uc_:\n"); uc_.print();
    // Rprintf("rr_:\n"); rr.print();
	
    K(u, u) = as_scalar(k_uu_upd2);
    // replace_uv_(K,  u__,  u__,   k_uu_upd2,  shift=shift); 

    K.submat(uc_, u_) = -K_uc_u_upd;
    // for (size_t k=0; k<uc_.n_elem; k++){
      // K(uc_(k),u) =  -K_uc_u_upd(uc_(k));
    // }
    // replace_uv_(K, -u__,  u__,  -K_uc_u_upd, shift=shift);

    K.submat(u_, uc_) = trans(-K_uc_u_upd);
    // replace_uv_(K,  u__, -u__,  -K_uc_u_upd, shift=shift);

    // as_colvec(RR3);    
    // replace_uv_(K, -u__, -u__,   RR3,        shift=shift);

    K.submat(uc_, uc_) = RR3;

// mat RR2 = CC2 + k_uu_upd2 * DD2 * DD2.t();
    // mat RR2 = CC2 + K_uc_u_upd * DD2.t();
    // Rprintf("RR2:\n"); RR2.print();    
    // Rprintf("K2 after smart update:\n"); K2.print();    
    // mat diff = K - K2;
    // Rprintf("diff:\n"); diff.print();
  } else {
    // Rprintf("brute force update of K\n");
    // solve_qr_(Sigma);
    K = inv_qr_(Sigma);
    // K = solve_qr_R(Sigma);    
  }
  // Rprintf("K after brute force update:\n"); K.print();
}


// [[Rcpp::export]]
mat update_row_K_(int u, mat& Sigma, mat& K, bool update_K=false, bool print=false){
  update_row_K2_(u, Sigma, K, update_K, print);
  return(K);
}


//[[Rcpp::export]]
mat update_row_Sigma_(int u, mat& Sigma, const mat& amat, bool print=false){
  update_row_Sigma2_(u, Sigma, amat, print);
  return(Sigma);
}


void innerloop1_update_Sigma2_(mat& Sigma, mat& amat){
  for (size_t u=0; u<amat.n_rows; u++){
    update_row_Sigma2_(u, Sigma, amat);
  }
}


void innerloop2_update_Sigma_K2_(mat& Sigma, mat& K, mat& amat, bool update_K=true, bool print=false){
  // Rprintf("innerloop2_update_Sigma_K2_\n");
  for (size_t u=0; u<amat.n_rows; u++){
    update_row_Sigma2_(u, Sigma, amat);
    update_row_K2_    (u, Sigma, K, update_K=update_K, print=print);
  }
}


//[[Rcpp::export]]
List outerloop1_(mat& Sigma, mat& K, umat& Emat, umat& Emat_c, mat& amat, int& nobs, double& eps, int& maxit){
  int it1 = 0;
  bool converged = false;
    
  double mad, conv_crit;
  // double logLp = ggm_logL_(Sigma, K, nobs);
  
  mat Sigma_prev = diagmat(Sigma.diag());
  // Rprintf("Sigma_prev:\n"); Sigma_prev.print();

  while (!converged){
    innerloop1_update_Sigma2_(Sigma, amat);
    mad = mean_abs_diff_non_edge_(Sigma, Sigma_prev, Emat_c); // FIXME for testing
    Sigma_prev = Sigma;
    conv_crit = mad;
    // Rprintf("conv_crit %f\n", conv_crit);
      
    it1++;
    if ((it1 == maxit) || (conv_crit < eps)){ break;}
  }

  int itcount= it1;
  return List::create(_["iter"]  = itcount, _["mad"]=mad); //FIXME mad should be conv_crit		
}



List outerloop2_(mat& Sigma, mat& K, umat& Emat, umat& Emat_c, mat& amat, double& eps, int& maxit){

  bool update_K = true;
  
  double dif2, conv_crit;

  double d = det(Sigma);

  if (d > 0){
    K = inv_qr_(Sigma);
    // K = eye(Sigma.n_rows, Sigma.n_cols);
  } else {
    stop("NCD algorithm failed");
  }

  // Rprintf("Sigma, K and KSigma (before updating)\n");
  // mat KSig = K * Sigma;  Sigma.print(); K.print(); KSig.print();

  // diff = KSig - eye(K.n_rows, K.n_cols);
  
  bool converged = false;
  int it2 = 0;
  while (!converged){
    innerloop2_update_Sigma_K2_(Sigma, K, amat, update_K=update_K);
    dif2 = diff_fun_(Sigma, K, Emat_c); // FIXME for testing
    it2++;
    conv_crit = dif2;
    if ((it2 == maxit) || (conv_crit < eps)){ break;}
  }

  int itcount= it2;
  return List::create(_["iter"] = itcount, _["conv_crit"] = conv_crit);		

}

//[[Rcpp::export(.c_ncd_ggm_)]]
List ncd_ggm_(mat& S, List& Elist, umat& Emat, int& nobs,
	       mat K,       
	       int& iter, double& eps, int& convcrit, int& print, List& aux){

  mat Sigma  = S;
  // Rprintf("emat:\n"); Emat.print();
  mat amat   = as_emat2amat_(Emat-1, S.n_rows);
  // Rprintf("amat:\n"); amat.print();
  umat Emat_c = as_emat_complement_(Emat - 1, S.n_rows);
  // Rprintf("emat_c:\n"); Emat_c.print();
  
  List res1 = outerloop1_(Sigma, K, Emat, Emat_c, amat, nobs, eps, iter);
  // Rprintf("Sigma after outerloop1:\n"); Sigma.print();

  // Rprintf("outerloop1 iter = %d\n", (int) res1["iter"]);
  
  List res2 = outerloop2_(Sigma, K, Emat, Emat_c, amat, eps, iter);
  // Rprintf("Sigma and K after outerloop2 :\n"); Sigma.print(); K.print();
  // Rprintf("outerloop2 iter = %d\n", (int) res2["iter"]);
  
  // Rprintf("Sigma:\n"); Sigma.print();
  // Rprintf("K:\n"); K.print();

  int itcount = (int) res2["iter"] + (int) res1["iter"];
  double conv_check = res2["conv_crit"];
  double logL = ggm_logL_(S, K, nobs);  

  RETURN_VALUE;
}



