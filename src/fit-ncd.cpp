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


double duality_gap_(mat& S, mat& K, int nobs){
  mat KS = K * S;
  double out = nobs * (accu(KS) - log(det(KS)) - S.n_rows) / 2; // FIXME Fragile
  return out;
}

// ### CALLING R FUNCTIONS


double mean_abs_diff_non_edge_(mat& Sigma1, mat& Sigma2, umat& emc){
  double out = accu(abs(Sigma1 - Sigma2)) / (2 * emc.n_cols);
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

void update_row_Sigma_(int u, mat& Sigma, const mat& amat, int nobs, int print=0){


  // FIXME Need not be computed each time...
  uvec u_    = {(unsigned int) u};      // convert int to uvec
  uvec ne_u_ = find(amat.rows(u_) > 0); // Returns column vector
  int deg_u  = accu(amat.rows(u_));

  if (print >= 4){
    Rprintf("++++ Updating Sigma for u=%i with degree %i\n", u, deg_u);
  }

  vec w(Sigma.n_cols, fill::zeros);  
  vec tt;
  
  if (ne_u_.n_rows > 0){  
    mat s12 = Sigma.cols(u_);
    mat AA = Sigma.submat(ne_u_, ne_u_);
    vec bb = s12.rows(ne_u_);

    if (deg_u > nobs){
      tt = pinv(AA) * bb;       // Rprintf("using pinv\n");
    } else {
      tt = solve(AA, bb);       // Rprintf("using solve\n");      
    }
      
    w = Sigma.cols(ne_u_) * tt;
  }

  // Rprintf("inserting:\n");
  double sigma_uu = Sigma(u, u); // Store this because element is overwritten below    
  Sigma.col(u) = w;
  Sigma.row(u) = w.t();
  Sigma(u, u)  = sigma_uu;       // Restore element
}


void update_row_K_(int u, mat& Sigma, mat& K, int smart=0, double eps_smart=0.0, int print=0){
  
  if (print >= 4){
    if (smart == 3)
      Rprintf("++++ Updating (brute force) K for u=%i\n", u);
    else
      Rprintf("++++ Updating (smart) K for u = %3i, smart=%i \n", u, smart);
  }

  if (smart==3){
    K = inv_qr_(Sigma);    
  } else {

    uvec u_  = {(unsigned int) u};    
    ivec u__ = {(int) u};

    arma::uvec uc_  = arma::linspace<arma::uvec>(0, K.n_cols - 1, K.n_cols);
    uc_.shed_rows(u_);

    double k_uu     = as_scalar(K(u, u));    
    double sigma_uu = as_scalar(Sigma(u, u));

    mat K_uc_u      = K.submat(uc_, u_);
    mat K_uc_uc     = K.submat(uc_, uc_);
    mat Sigma_uc_u  = Sigma.submat(uc_, u_);
    
    mat CC2 = K_uc_uc - K_uc_u * (trans(K_uc_u) / k_uu);
    mat DD2 = CC2 * Sigma_uc_u;

    // DO UPDATE
    mat k_uu_upd2  = 1 / (sigma_uu - trans(Sigma_uc_u) * DD2);
    mat K_uc_u_upd = as_scalar(k_uu_upd2) * DD2;
    	
    K(u, u) = as_scalar(k_uu_upd2);
    K.submat(uc_, u_) = -K_uc_u_upd;
    K.submat(u_, uc_) = trans(-K_uc_u_upd);

    mat new2, old2, RR3, dd;
    double ee;
      
    if ((smart == 1) || (smart == 2)){
      new2 = K_uc_u_upd * (trans(K_uc_u_upd) / as_scalar(k_uu_upd2));
      old2 = K_uc_u     * (trans(K_uc_u)     / as_scalar(k_uu));
      // Rprintf("new2:\n"); new2.print();
      // Rprintf("old2:\n"); old2.print();
    }

    switch (smart) {
    case 0: // The glasso type update
      break;
    case 1: // Update as described in paper
      RR3  = K_uc_uc + new2 - old2;
      K.submat(uc_, uc_) = RR3;
      break;
    case 2: // Approximate updates
      dd = K_uc_u_upd - K_uc_u;
      // Rprintf("dd:\n"); dd.print();
      ee = accu(dd.t() * dd) / Sigma.n_rows;      
      // Rprintf("ee: %d\n", ee);
      if (ee >= eps_smart){
    	if (print >= 4)
	  Rprintf("++++ Updating u=%3i ee = %f\n", u, ee);
    	RR3  = K_uc_uc + new2 - old2;
    	K.submat(uc_, uc_) = RR3;
      }
      break;
    default:
      Rcout << "Invalid value of smart\n";
    }    
  } 
}



void innerloop2_update_Sigma_K_(mat& Sigma, mat& K, mat& amat, int nobs,
				int smart=0, double eps_smart=0.0, int print=0){
  if (print >= 3){
    Rprintf("+++ Running innerloop2_update_Sigma_K\n");
   }

  for (size_t u=0; u<amat.n_rows; u++){
    update_row_Sigma_(u, Sigma, amat, nobs, print);
    update_row_K_    (u, Sigma, K, smart=smart, eps_smart=eps_smart, print=print);
  }
}

void innerloop1_update_Sigma_(mat& Sigma, mat& amat, int nobs, int print=0){

  if (print >= 3){
    Rprintf("+++ Running innerloop1_update_Sigma\n");
   }
  
  for (size_t u=0; u<amat.n_rows; u++){
    update_row_Sigma_(u, Sigma, amat, nobs, print);
  }
}



List outerloop2_(mat& Sigma, mat& K, umat& Emat, umat& Emat_c, mat& amat, double& eps, int& maxit,
		 int rank_Sigma,
		 int smart=0, double eps_smart=0.0, int print=0){

  if (print >=2){
    if (smart==3) Rprintf("++ Running outerloop2 - brute force update of K\n");
    else          Rprintf("++ Running outerloop2 - smart update of K, smart=%d, eps_smart=%f\n", smart, eps_smart);
  }

  double dif2, conv_crit;
  bool converged = false;
  int it2 = 0;
  while (!converged){
    innerloop2_update_Sigma_K_(Sigma, K, amat, smart=smart, eps_smart=eps_smart, print=print);
    dif2 = diff_fun_(Sigma, K, Emat_c); // FIXME for testing
    it2++;
    conv_crit = dif2;
    if (smart == 0) {
      // Rprintf("smart=0; exiting..\n");
      break;}
    if ((it2 == maxit) || (conv_crit < eps)){ break;}
  }
  return List::create(_["iter"] = it2, _["conv_crit"] = conv_crit);		
}



//[[Rcpp::export]]
List outerloop1_(mat& Sigma, mat& K, umat& Emat, umat& Emat_c, mat& amat,
		 int& nobs, double& eps, int& maxit, int print=0){
  int it1 = 0;
  bool converged = false;

  if (print >=2){
    Rprintf("++ Running outerloop1\n");
  }
    
  double mad, conv_crit;
  mat Sigma_prev = diagmat(Sigma.diag());
  
  // double logLp = ggm_logL_(Sigma, K, nobs);
  // Rprintf("Sigma_prev:\n"); Sigma_prev.print();

  while (!converged){
    innerloop1_update_Sigma_(Sigma, amat, nobs, print);
    mad = mean_abs_diff_non_edge_(Sigma, Sigma_prev, Emat_c); // FIXME for testing

    Sigma_prev = Sigma;
    conv_crit = mad;

    it1++;
    if ((it1 == maxit) || (conv_crit < eps)){ break;}
  }

  return List::create(_["iter"]  = it1, _["mad"]=mad); //FIXME mad should be conv_crit		
}



//[[Rcpp::export(.c_ncd_ggm_)]]
List ncd_ggm_(mat& S, List& Elist, umat& Emat, int& nobs,
	       mat K,       
	       int& iter, double& eps, int& convcrit, int print, List& aux){

  int  smart       = aux["smart"];
  double eps_smart = aux["eps_smart"];
  
  umat Emat_c = as_emat_complement_(Emat - 1, S.n_rows);
  mat amat    = as_emat2amat_(Emat-1, S.n_rows);
  mat Sigma   = S;


  List res1 = outerloop1_(Sigma, K, Emat, Emat_c, amat, nobs, eps, iter, print);

  uword rank_Sigma = arma::rank(Sigma);  
  if (rank_Sigma >= Sigma.n_cols){
    K = inv_qr_(Sigma);
  } else {
    REprintf("Rank of Sigma = %d nobs = %d\n", rank_Sigma, nobs);
    stop("NCD algorithm not convergent\n");
  }
  
  List res2 = outerloop2_(Sigma, K, Emat, Emat_c, amat, eps, iter,
			  rank_Sigma=rank_Sigma,
			  smart=smart, eps_smart=eps_smart, print=print);

  double gap = duality_gap_(Sigma, K, nobs);
  int itcount = (int) res2["iter"] + (int) res1["iter"];
  double conv_check = res2["conv_crit"];
  double logL = ggm_logL_(S, K, nobs);  

  RETURN_VALUE;
}






  // if (print >= 5){
  //   Rprintf("+++++ emat:\n"); Emat.print();
  //   Rprintf("+++++ amat:\n"); amat.print();
  //   Rprintf("+++++ emat_c:\n"); Emat_c.print();
  // }




































  // Rprintf("Sigma after outerloop1:\n"); Sigma.print();
  // Rprintf("outerloop1 iter = %d\n", (int) res1["iter"]);


  // Rprintf("Sigma and K after outerloop2 :\n"); Sigma.print(); K.print();
  // Rprintf("outerloop2 iter = %d\n", (int) res2["iter"]);
  
  // Rprintf("Sigma:\n"); Sigma.print();
  // Rprintf("K:\n"); K.print();


    // replace_uv_(K,  u__,  u__,   k_uu_upd2,  shift=shift); 
    // replace_uv_(K,  u__, -u__,  -K_uc_u_upd, shift=shift);
    
    // Rprintf("K_uc_u:\n"); K_uc_u.t().print();
    // Rprintf("K_uc_uc:\n"); K_uc_uc.print();
    // Rprintf("Sigma_uc_u:\n"); Sigma_uc_u.t().print();

    // Rprintf("CC2:\n"); CC2.print();
    // Rprintf("DD2:\n"); DD2.t().print(); 

    // for (size_t k=0; k<uc_.n_elem; k++){
      // K(uc_(k),u) =  -K_uc_u_upd(uc_(k));
    // }
    // replace_uv_(K, -u__,  u__,  -K_uc_u_upd, shift=shift);


    // arma::uvec vv_  = arma::linspace<arma::uvec>(0, K.n_cols - 1, K.n_cols);
    // arma::uvec uc_  = setdiff_(vv_, u_);
    // arma::uvec rr = vv_;
    // rr.shed_rows(u_);
    // Rprintf("uc_:\n"); uc_.print();
    // Rprintf("rr_:\n"); rr.print();


    // mat K_uc_u     = extract_uv_(K, -u__,  u__, shift=shift);
    // mat K_uc_uc    = extract_uv_(K, -u__, -u__, shift=shift);
    // mat Sigma_uc_u = extract_uv_(Sigma, -u__, u__, shift=shift);


    
// mat RR2 = CC2 + k_uu_upd2 * DD2 * DD2.t();
    // mat RR2 = CC2 + K_uc_u_upd * DD2.t();
    // Rprintf("RR2:\n"); RR2.print();    
    // Rprintf("K2 after smart update:\n"); K2.print();    
    // mat diff = K - K2;
    // Rprintf("diff:\n"); diff.print();



      // NOT A GOOD IDEA TIMEWISE
      // double k_uu_upd2_inv = 1 / as_scalar(k_uu_upd2);
      // double k_uu_inv      = 1 / as_scalar(k_uu);
      // new2 = (k_uu_upd2_inv * K_uc_u_upd) * trans(K_uc_u_upd);
      // old2 = (k_uu_inv      * K_uc_u)     * trans(K_uc_u);



    // Rprintf("AA:\n"); AA.print();
    // Rprintf("bb:\n"); bb.t().print();
    // Rprintf("tt:\n"); tt.t().print();
    
    // vec beta(Sigma.n_cols, fill::zeros);
    // beta.elem(ne_u_) = tt;
    // w = Sigma * beta;
    // Rprintf("beta:\n"); beta.t().print();
