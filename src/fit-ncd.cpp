#include "RcppArmadillo.h"
#include "grips-utils.h"
#include "general-utils.h"
#include "arma_utils.h"
#include "convergence.h"
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


#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)>(b))?(b):(a))


// ### ###################################################
// ### UTILITY FUNCTIONS FOR NCD
// ### ###################################################

double duality_gap_(mat& S, mat& K, int nobs){
  mat KS = K * S;

  double val, sign;
  log_det(val, sign, KS);

  // double out = nobs * (accu(K % S) - log(det(KS)) - S.n_rows) / 2; // FIXME Fragile
  double out = nobs * (accu(K % S) - val - S.n_rows) / 2; // FIXME Fragile
  return out;
}

double mean_abs_diff_non_edge_(mat& Sigma1, mat& Sigma2, umat& emc){
  double out = accu(abs(Sigma1 - Sigma2)) / (2 * emc.n_cols);
  return out;
}

double max_abs_diff_non_edge_(mat& Sigma1, mat& Sigma2, umat& emc){
  double out = max(max(abs(Sigma1 - Sigma2)));
  return out;
}

// [[Rcpp::export]]
double diff_fun_(mat& Sigma, mat& K, umat emc){
  uvec ind = sub2ind(size(Sigma), emc);
  vec  Kuv = K(ind);
  vec  diagS = Sigma.diag();
  uvec r0 = {0}, r1 = {1};
  mat  emc2 = conv_to<mat>::from(emc);
  uvec s0 = conv_to<uvec>::from(emc2.rows(r0));
  uvec s1 = conv_to<uvec>::from(emc2.rows(r1));
  // Rprintf("s0\n");  s0.print();
  // Rprintf("s1\n");  s1.print();
  vec Suu = diagS.elem(s0);
  vec Svv = diagS.elem(s1);
  // Rprintf("Suu\n");  Suu.print();
  // Rprintf("Svv\n");  Svv.print();
  // int d = Sigma.n_rows * (Sigma.n_rows + 1) / 2;
  // double out = sum(abs(abs(Kuv) % sqrt(Suu % Svv))) / d;
  double out = max(max(abs(Kuv) % sqrt(Suu % Svv)));  
  return out;
}



// ### ###################################################
// ### NCD FUNCTIONS
// ### ###################################################



// ### ###################################################
// ### update rows of Sigma and K
// ### ###################################################

// *** Used by outerloop1 and outerloop2

void update_Sigma_row_(int u, mat& Sigma, const mat& amat, int nobs, int print=0){

  // FIXME Need not be computed each time...
  uvec u_    = {(unsigned int) u};      // convert int to uvec
  uvec ub_   = find(amat.rows(u_) > 0); // Returns column vector
  int  deg_u = accu(amat.rows(u_));

  if (print >= 4){
    Rprintf("++++ Updating Sigma for u=%i with degree %i\n", u, deg_u);
  }

  vec w(Sigma.n_cols, fill::zeros);  
  vec beta_star;
  
  if (ub_.n_rows > 0){  
    mat AA    = Sigma.submat(ub_, ub_);
    mat s_u   = Sigma.cols(u_);
    vec s_ubu = s_u.rows(ub_);

    if (deg_u > nobs - 1){
      beta_star = pinv(AA) * s_ubu;       // Rprintf("using pinv\n");
    } else {
      beta_star = solve(AA, s_ubu);       // Rprintf("using solve\n");      
    }
    w = Sigma.cols(ub_) * beta_star;
  }

  // Rprintf("inserting:\n");
  double sigma_uu = Sigma(u, u); // Store this because element is overwritten below    
  Sigma.col(u) = w;
  Sigma.row(u) = w.t();
  Sigma(u, u)  = sigma_uu;       // Restore element
}


// *** Used by outerloop2

void update_K_row_(int u, mat& Sigma, mat& K, const mat& amat, int smart=0, double eps_smart=0.0, int print=0){
  
  if (smart==3){
    K = inv_qr_(Sigma);    
  } else {
    uvec u_    = {(unsigned int) u};    
    uvec ub_   = find(amat.rows(u_) > 0); // Returns column vector
    uvec uc_   = arma::linspace<arma::uvec>(0, K.n_cols - 1, K.n_cols);
    int  deg_u = accu(amat.rows(u_));
    uc_.shed_rows(u_);

    // ivec u__ = {(int) u};
    if (print >= 4){
      Rprintf("++++ Updating K for u=%i with degree %i\n", u, deg_u);
      Rprintf("++++ ub_: ");  ub_.t().print();
    }

    double k_uu     = as_scalar(K(u, u));    
    double sigma_uu = as_scalar(Sigma(u, u));

    mat K_ucu      = K.submat(uc_, u_);
    mat K_ucuc     = K.submat(uc_, uc_);
    mat Sigma_ucu  = Sigma.submat(uc_, u_);
    mat CC2        = K_ucuc - K_ucu * (trans(K_ucu) / k_uu);
    mat DD2        = CC2 * Sigma_ucu;
    mat k2_uu      = 1 / (sigma_uu - trans(Sigma_ucu) * DD2);
    mat K2_ucu     = as_scalar(k2_uu) * DD2;

    // DO UPDATE    
    mat RR2, RR, K2_ucuc;

    K(u, u)           = as_scalar(k2_uu);
    K.submat(uc_, u_) = -K2_ucu;
    K.submat(u_, uc_) = trans(-K2_ucu);
      
    RR2 = K2_ucu * (trans(K2_ucu) / as_scalar(k2_uu));
    RR  = K_ucu  * (trans(K_ucu)  / as_scalar(k_uu));
    K2_ucuc  = K_ucuc + RR2 - RR;
    K.submat(uc_, uc_) = K2_ucuc;


    // mat K_ubu      = K.submat(ub_, u_);
    // mat K_ubub     = K.submat(ub_, ub_);
    // mat Sigma_ubu  = Sigma.submat(ub_, u_);
    // mat CC3        = K_ubub - K_ubu * (trans(K_ubu) / k_uu);
    // mat DD3        = CC3 * Sigma_ubu;
    // mat k3_uu      = 1 / (sigma_uu - trans(Sigma_ubu) * DD3);
    // mat K3_ubu     = as_scalar(k3_uu) * DD3;

    // mat AA3, AA, K3_ubub;

    // K(u, u)           = as_scalar(k3_uu);
    // K.submat(ub_, u_) = -K3_ubu;
    // K.submat(u_, ub_) = trans(-K3_ubu);

    // AA3 = K3_ubu * (trans(K3_ubu) / as_scalar(k3_uu));
    // AA  = K_ubu  * (trans(K_ubu)  / as_scalar(k_uu));  
    // K3_ubub  = K_ubub + AA3 - AA;
    // K.submat(ub_, ub_) = K3_ubub;        

  } 
}




// ### ###################################################
// ### OUTERLOOP1
// ### ###################################################


void innerloop1_update_Sigma_(mat& Sigma, mat& amat, int nobs, int print=0){

  if (print >= 4){
    Rprintf("++++ Running innerloop1_update_Sigma\n");
  }
  
  for (size_t u=0; u<amat.n_rows; u++){
    update_Sigma_row_(u=u, Sigma=Sigma, amat=amat, nobs=nobs, print=print);
  }
}

//[[Rcpp::export]]
List outerloop1_(mat& Sigma, mat& K, umat& emat, umat& emat_c, mat& amat,
		 int& nobs, double& eps, int maxiter, int print=0){

  if (print >=2){
    Rprintf("++ Running outerloop1\n");
  }
  
  bool converged = false;    
  double mad, conv_crit, mno;
  mat Sigma_prev = diagmat(Sigma.diag());
  int iter = 0;

  // Rcout << "eps :" << eps << " n_rows : " << Sigma.n_rows << " eps2 : " << eps2 <<"\n" ;
  while (!converged){
    innerloop1_update_Sigma_(Sigma=Sigma, amat=amat, nobs=nobs, print=print);
    // mad = mean_abs_diff_non_edge_(Sigma, Sigma_prev, emat_c); // FIXME for testing
    // mad = max_abs_diff_non_edge_(Sigma, Sigma_prev, emat_c); // FIXME for testing
    mat Delta = Sigma - Sigma_prev;
    mno = mnormone_(Delta);
    // Rprintf("iter: %d eps2: %14.10f mno : %14.10f\n", iter, eps2, mno);
    Sigma_prev = Sigma;
    conv_crit = mno;
    iter++;
    if ((iter == maxiter) || (conv_crit < eps)){
      break;
    }
  }

  return List::create(_["iter"]  = iter, _["mad"]=mno); //FIXME mad should be conv_crit		
}


// ### ###################################################
// ### OUTERLOOP2
// ### ###################################################

void innerloop2_update_Sigma_K_(mat& Sigma, mat& K, mat& amat, int nobs,
				int smart=0, double eps_smart=0.0, int print=0){
  if (print >= 3){
    Rprintf("+++ Running innerloop2_update_Sigma_K\n");
   }

  for (size_t u=0; u<amat.n_rows; u++){
    update_Sigma_row_(u=u, Sigma=Sigma,      amat=amat, nobs=nobs, print=print);    
    update_K_row_    (u=u, Sigma=Sigma, K=K, amat=amat, smart=smart, eps_smart=eps_smart, print=print);
  }
}

List outerloop2_(mat& Sigma, mat& K, umat& emat, umat& emat_c, mat& amat, int nobs, double& eps, int& maxiter,
		 int rank_Sigma,
		 int smart=0, double eps_smart=0.0, int print=0){

  if (print >=2){
    if (smart==3) Rprintf("++ Running outerloop2 - brute force update of K\n");
    else          Rprintf("++ Running outerloop2 - smart update of K, smart=%d, eps_smart=%f\n", smart, eps_smart);
  }

  double dif2, mno, conv_crit;
  bool converged = false;
  int iter = 0;
  while (!converged){
    innerloop2_update_Sigma_K_(Sigma=Sigma, K=K, amat=amat, nobs=nobs,
			       smart=smart, eps_smart=eps_smart, print=print);
    
    // dif2 = diff_fun_(Sigma, K, emat_c); // FIXME for testing
    // conv_crit = dif2;

    mat Delta = K - project_K_onto_G_(K, emat_c);
    mno = mnormone_(Delta);
    conv_crit = mno;

    if (print>=3)
      Rprintf("+++ innerloop2 iter : %d mno : %f\n", iter, mno);
    
    iter++;
    if ((iter == maxiter) || (conv_crit < eps)){
      break;
    }
  }
  return List::create(_["iter"] = iter, _["conv_crit"] = conv_crit);		
}


// ### ###################################################
// ### From K to Sigma
// ### ###################################################

void Sigma_to_K_row_(int u, mat& Sigma, mat& K, const mat& amat, int nobs, int print=0){

  // FIXME Need not be computed each time...
  uvec u_    = {(unsigned int) u};      // convert int to uvec
  uvec ub_   = find(amat.rows(u_) > 0); // Returns column vector
  int  deg_u = accu(amat.rows(u_));

  if (print >= 4){
    Rprintf("++++ Updating Sigma for u=%i with degree %i\n", u, deg_u);
  }
  
  vec beta_pad(Sigma.n_cols, fill::zeros);  
  vec beta_star, K_ubu;
  mat s_u, s_ubu, AA;
  
  if (ub_.n_rows > 0){  
    AA    = Sigma.submat(ub_, ub_);
    s_u   = Sigma.cols(u_);
    s_ubu = s_u.rows(ub_);
    
    if (deg_u > nobs - 1){
      beta_star = pinv(AA) * s_ubu;       // Rprintf("using pinv\n");
    } else {
      beta_star = solve(AA, s_ubu);       // Rprintf("using solve\n");      
    }
    beta_pad.elem(ub_) = beta_star;
  }
  
  double k2_uu = 1 / as_scalar(Sigma(u, u) - accu(s_ubu % beta_star));
  K_ubu = - beta_star * k2_uu;
  
  K.submat(ub_, u_) = K_ubu;
  K.submat(u_, ub_) = trans(K_ubu);
  K(u, u) = k2_uu;
}


//[[Rcpp::export]]
void Sigma_to_K_(mat& Sigma, mat& K, mat& amat, int nobs, int print=0){

  for (size_t u=0; u<amat.n_rows; u++){
    Sigma_to_K_row_(u=u, Sigma=Sigma, K=K, amat=amat, nobs=nobs, print=print);    
  }
}


// ### ###################################################
// ### MAIN NCD FUNCTION 
// ### ###################################################

#define CHECK_K								\
  uword rank_Sigma = arma::rank(Sigma, sqrt(datum::eps) * S.n_cols);	\
  if (rank_Sigma >= Sigma.n_cols){					\
    K = inv_qr_(Sigma);							\
  } else {								\
    REprintf("Rank of Sigma = %d nobs = %d\n", rank_Sigma, nobs);	\
    stop("NCD algorithm not convergent\n");				\
  }									\

bool has_full_rank(mat& Sigma){
  uword rank_Sigma = arma::rank(Sigma, sqrt(datum::eps) * Sigma.n_cols);
  return (rank_Sigma >= Sigma.n_cols); 	    
}

bool is_pos_def(mat& Sigma){
  uword rank_Sigma = arma::rank(Sigma, sqrt(datum::eps) * Sigma.n_cols);
  return (rank_Sigma >= Sigma.n_cols); 	    
}

//[[Rcpp::export(.c_ncd_ggm_)]]
List ncd_ggm_(mat& S, List& elist, umat& emat, int& nobs,
	      mat K,       
	      int maxiter, double& eps, int& convcrit, int print, List& aux){
  
  int version      = aux["version"];
  int smart        = aux["smart"];
  double eps_smart = aux["eps_smart"];
  
  umat emat_c = as_emat_complement_(emat-1, S.n_rows);
  mat amat    = as_emat2amat_(emat-1, S.n_rows);
  mat Sigma   = S, K2, Delta;
  List res1, res2;
  double logL, gap, conv_check;
  int iter1, iter2, itcount;

  double eps2 = MIN(eps, 1.0/Sigma.n_rows);  
  double mno;

  // FIXME NOTICE SCALING OG EPS2
  res1 = outerloop1_(Sigma=Sigma, K=K, emat=emat, emat_c=emat_c, amat=amat,
		     nobs=nobs, eps=eps2/100, maxiter=maxiter, print=print);
  
  iter1 = res1["iter"];
  if (print>=2)
    Rprintf("++ outerloop1 iterations : %d\n", iter1);
  
  CHECK_K;
  
  switch (version){    
  case 1: 
    // Rprintf("version 1 - full\n");
    
    // Check that Sigma has full rank
    if (!has_full_rank(Sigma)){
      // logL = NA; K=NA, dgap=NA
      // abort
      ;
    } else {    
      K = inv_qr_(Sigma);
      res2 = outerloop2_(Sigma=Sigma, K=K, emat=emat, emat_c=emat_c, amat=amat, nobs=nobs, eps=eps2, maxiter=maxiter,
			 rank_Sigma=rank_Sigma,
			 smart=smart, eps_smart=eps_smart, print=print);
      iter2 = res2["iter"];
      if (print>=2)
	Rprintf("++ outerloop2 iterations : %d\n", iter2);
      
      K2 = project_K_onto_G_(K, emat_c);
      Delta = K - K2;
      mno = mnormone_(Delta);
      if (print>=3)
	Rprintf("+++ full mno : %f\n", mno);
      conv_check = mno;
      if (iter2 < maxiter){ // Then K is posdef	
      	logL = ggm_logL_(S, K2, nobs);
      	gap  = duality_gap_(Sigma, K2, nobs);
      } else {
      	REprintf("Algorithn may not have converged\n");
      	// K = NA; upper_limit_logL = formel (23)
      }
      itcount = iter1 + iter2;
    }
    break;
      
  case 0:
    // Rprintf("version 0 - fast\n");
    if (!has_full_rank(Sigma)){
      // logL = NA; K=NA, dgap=NA
      // abort
    } else {    
      K = inv(Sigma);
      K2 = project_K_onto_G_(K, emat_c);
      Delta = K - K2;
      mno = mnormone_(Delta);
      if (print>=3)
	Rprintf("+++ fast mno : %f\n", mno);
      conv_check = mno;      
      if (!is_pos_def(K2)){
	REprintf("Algorithm may not have converged\n");
	// K = NA; upper_limit_logL = formel (23)
      } else {
	logL = ggm_logL_(S, K2, nobs);
	gap  = duality_gap_(Sigma, K2, nobs);	
      }
      itcount = iter1 + 1;
    }
    break;
    
  default:
    Rprintf("'version' must be 0, 1\n");
  }

  RETURN_VALUE;
}




























// //[[Rcpp::export(.c_ncd_ggm_)]]
// List ncd_ggm_(mat& S, List& Elist, umat& emat, int& nobs,
// 	       mat K,       
// 	       int iter, double& eps, int& convcrit, int print, List& aux){

//   int version      = aux["version"];
//   int smart        = aux["smart"];
//   double eps_smart = aux["eps_smart"];
  
//   umat emat_c = as_emat_complement_(emat - 1, S.n_rows);
//   mat amat    = as_emat2amat_(emat-1, S.n_rows);
//   mat Sigma   = S;
//   List res1, res2;
//   double gap, conv_check;
//   int maxit, itcount;
  
//   res1 = outerloop1_(Sigma=Sigma, K=K, emat=emat, emat_c=emat_c, amat=amat,
// 		     nobs=nobs, eps=eps, maxit=iter, print=print);  
//   CHECK_K;
    
//   switch (version){

//   case 1: // The glasso type update
//     Rprintf("version 1\n");
//     conv_check = res1["mad"];    
//     K = mat(S.n_rows, S.n_rows, fill::eye);
//     Sigma_to_K_(Sigma=Sigma, K=K, amat=amat, nobs=nobs, print=print);
//     // FIXME: Need sanity check of K being p.d.
//     gap = duality_gap_(Sigma, K, nobs);
//     itcount = (int) res1["iter"];    
//     break;

//   case 2:
//     Rprintf("version 2\n");
//     smart = 1; // Full update 
//     res2 = outerloop2_(Sigma=Sigma, K=K, emat=emat, emat_c=emat_c, amat=amat, nobs=nobs, eps=eps, maxit=iter,
// 		       rank_Sigma=rank_Sigma,
// 		       smart=smart, eps_smart=eps_smart, print=print);
//     K = mat(S.n_rows, S.n_rows, fill::eye);
//     Sigma_to_K_(Sigma=Sigma, K=K, amat=amat, nobs=nobs, print=print);
//     // FIXME: Need sanity check of K being p.d.
//     conv_check = res2["conv_crit"];
    
//     gap = duality_gap_(Sigma, K, nobs);
//     itcount = (int) res2["iter"] + (int) res1["iter"];
//     break;


//   case 3:
//     Rprintf("version 3\n");
//     smart = 1; // Full update 
//     res2 = outerloop2_(Sigma=Sigma, K=K, emat=emat, emat_c=emat_c, amat=amat, nobs=nobs, eps=eps, maxit=iter,
// 		       rank_Sigma=rank_Sigma,
// 		       smart=smart, eps_smart=eps_smart, print=print);    
//     conv_check = res2["conv_crit"];

//     K = project_K_onto_G_(K, emat_c);
//     // FIXME: Need sanity check of K being p.d.
//     gap = duality_gap_(Sigma, K, nobs);
//     itcount = (int) res2["iter"] + (int) res1["iter"];
//     break;
    
//   default:
//     Rprintf("'version' must be 1, 2 or 3\n");
//   }
 
//   double logL = ggm_logL_(S, K, nobs);  

//   RETURN_VALUE;
// }









    // RR.print();
    // AA.print();
    // Rprintf("K2_ucu:\n"); K2_ucu.t().print();
    // Rprintf("K_ucu:\n"); K_ucu.t().print();
    
    // switch (smart) {
    // case 0: // The kalif type update
    //   break;
    // case 1: // Update as described in paper
      
    //   break;
    // case 2: // Approximate updates
    //   Kdiff_ucu = K2_ucu - K_ucu;
    //   Kdiff_norm = accu(Kdiff_ucu.t() * Kdiff_ucu) / Sigma.n_rows;      
    //   if (Kdiff_norm >= eps_smart){
    // 	if (print >= 4)
    // 	  Rprintf("++++ Updating u=%3i Kdiff_norm = %f\n", u, Kdiff_norm);
    // 	K2_ucuc  = K_ucuc + RR2 - RR;
    // 	K.submat(uc_, uc_) = K2_ucuc;
    //   }
    //   break;
    // default:
    //   Rcout << "Invalid value of smart\n";
    // }


  // if (print >= 5){
  //   Rprintf("+++++ emat:\n"); emat.print();
  //   Rprintf("+++++ amat:\n"); amat.print();
  //   Rprintf("+++++ emat_c:\n"); emat_c.print();
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
    // beta.elem(ub_) = tt;
    // w = Sigma * beta;
    // Rprintf("beta:\n"); beta.t().print();
