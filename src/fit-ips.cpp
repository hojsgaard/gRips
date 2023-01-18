#include "RcppArmadillo.h"
#include "grips-utils.h"
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



// -------------------------------------------------------------
// ------ IPS - Classical iterative proportional scaling -------
// -------------------------------------------------------------

// Create list with complements to components in edges.
// [[Rcpp::export]]
List make_clist_(arma::mat& S, List& edges){
  int_vec idx = seq(1, S.n_rows);
  // Rcout << idx << std::endl;
  List clist = List(edges.length());
  for (int i=0; i<edges.length(); i++){
    int_vec cc = edges[i];
    int_vec aa = setdiff(idx, cc);
    clist[i] = aa;
  }
  // Rf_PrintValue(wrap(edges));   Rf_PrintValue(wrap(clist));
  return clist;
}

//[[Rcpp::export]]
List Scc_inv_list_(const mat& S, const List& edges0){
  List out(edges0.length());
  for (int i=0; i<edges0.length(); i++){
    uvec cc = as<arma::uvec>(edges0[i]);
    mat Scc = S.submat(cc, cc); 
    out[i] = inv(Scc);
  }
  return out;
}

//[[Rcpp::export]]
List Scc_list_(const mat& S, const List& edges0){
  List out(edges0.length());
  for (int i=0; i<edges0.length(); i++){
    uvec cc = as<arma::uvec>(edges0[i]);
    mat Scc = S.submat(cc, cc); 
    out[i] = Scc;
  }
  return out;
}

inline void ips_ggm_update_cc_parm_(mat& S, uvec& cc0, mat& K, uvec& aa0, mat& Scc_inv)
{
  // mat Scc=S.submat(cc0, cc0);
  mat Kaa=K.submat(aa0, aa0);
  mat Kca=K.submat(cc0, aa0);    
  mat Kac=K.submat(aa0, cc0);    
  // mat L = Kca * inv(Kaa) * Kac;
  K.submat(cc0, cc0) = Scc_inv + Kca * inv(Kaa) * Kac;
}

void ips_inner_(arma::mat& S, List& edges0,
		arma::mat& K,
		List& clist0){  // NOTE: edges0, clist0 are 0-based

  List Scc_inv_list = Scc_inv_list_(S, edges0);

  // Rcout << "ips_inner started" << std::endl;
  for (int i=0; i < edges0.length(); i++)
    {    
      uvec cc0 = (edges0[i]), aa0 = (clist0[i]);
      mat Scc_inv = Scc_inv_list[i];
      ips_ggm_update_cc_parm_(S, cc0, K, aa0, Scc_inv);
    } 
  // Rcout << "ips_inner done" << std::endl;
}


//[[Rcpp::export(.c_ips_ggm_)]]
List ips_ggm_(arma::mat& S, List& Elist, umat& Emat, int& nobs,
	      arma::mat K,
	      int& iter, double& eps, int& convcrit, int& print, List& aux){

  double logL, logLp, mad, conv_check=9999, conv_ref;
  char buffer[200];
  int itcount  = 0;
  double nparm = S.n_cols + Emat.n_cols;  
  umat Emat0   = Emat - 1;

  // Variables for ips
  mat Sigma = inv(K);
  List Elist0 = clone(Elist);
  List clist0 = make_clist_(S, Elist);
  
  for (int i=0; i<Elist.length(); i++){
    Elist0[i] = as<arma::uvec>(Elist0[i]) - 1; // 0-based
    clist0[i] = as<arma::uvec>(clist0[i]) - 1; // 0-based			     
  }
  // END
  
  INIT_CONVERGENCE_CHECK;
    
  for (; itcount < iter; ){  
    ips_inner_(S, Elist0, K, clist0);
    ++itcount;          

    if (true){ 
      switch(convcrit){
      case 1: 
	Sigma = inv(K); // Needed for ips only
	mad   = mean_abs_diff_on_Emat_(S, Sigma, Emat0, 0);
	conv_check = mad;
	PRINT_CONV_CHECK1;
	break;
      case 2:
	CONV_CHECK_LOGL_DIFF;
	break;
      case 3:
	CONV_CHECK_LOGL_DIFF_REF;	
	PRINT_CONV_CHECK3;
	break;
      }    
      if (conv_check < eps) break;
    }
  }

  logL  = ips_logL_(S, K, nobs);
  Sigma = inv(K);
  RETURN_VALUE;
}

  

// ------------------------------------------------------------------
// ------ FIPS - Fast iterative proportional scaling          -------
// ------------------------------------------------------------------


void fips_ggm_update_cc_parm0_(const mat& Scc, const uvec& cc0, mat& K, mat& Sigma,
				      const mat& Scc_inv)
{
  
  mat Kcc      = K.submat(cc0, cc0);	
  mat Sigmacc  = Sigma.submat(cc0, cc0);
  
  mat Kstar  = inv(Sigmacc);
  mat Kupd   = Scc_inv + Kcc - Kstar;
  mat Haux   = Kstar - Kstar * Scc * Kstar;  // Was Saux

  K.submat(cc0, cc0) = Kupd;

  Sigma -= Sigma.cols(cc0) * Haux * Sigma.rows(cc0);
  // mat V = Sigma.cols(cc0);
  // Sigma -= V * Haux * V.t();
}


void fips_ggm_update_cc_parm_(const mat& Scc, const uvec& cc0, mat& K, mat& Sigma,
				      const mat& Scc_inv)
{
  
  mat Kcc      = K.submat(cc0, cc0);	
  mat Sigmacc  = Sigma.submat(cc0, cc0);
  
  mat Kstar  = inv(Sigmacc);
  
  mat Kupd   = Scc_inv + Kcc - Kstar;
  K.submat(cc0, cc0) = Kupd;
  // K.submat(cc0, cc0) = Scc_inv + Kcc - Kstar;

  mat Haux   = Kstar - Kstar * Scc * Kstar;  // Was Saux
  Sigma -= Sigma.cols(cc0) * Haux * Sigma.rows(cc0);
  // Sigma -= Sigma.cols(cc0) * (Kstar - Kstar * Scc * Kstar) * Sigma.rows(cc0);
}


void fips_inner_(const mat& S, const List& Elist0,
		  mat& K, mat& Sigma, List& Scc_list, List& Scc_inv_list)
{
  for (int i=0; i < Elist0.length(); ++i){
    uvec cc0 = Elist0[i];
    fips_ggm_update_cc_parm_(Scc_list[i], cc0, K, Sigma, Scc_inv_list[i]);
  }
}

//[[Rcpp::export(.c_fips_ggm_)]] 
List fips_ggm_(mat& S, List& Elist, umat& Emat, int& nobs,
	       mat K,       
	       int& iter, double& eps, int& convcrit, int& print, List& aux){

  double logL, logLp, mad, conv_check=9999, conv_ref;
  char buffer[200];
  int itcount  = 0;
  double nparm = S.n_cols + Emat.n_cols;
  umat Emat0   = Emat - 1;
 
  // Variables for fips
  mat Sigma    = initSigma_(S);  
  // END
  
  INIT_CONVERGENCE_CHECK;

  // Compute Scc etc only once
  List Elist0 = clone(Elist);  
  for (int i=0; i < Elist.length(); i++) {
    Elist0[i] = as<arma::uvec>(Elist0[i]) - 1; // 0-based
  }
  
  List Scc_inv_list = Scc_inv_list_(S, Elist0);
  List Scc_list     = Scc_list_(S, Elist0);

  for (; itcount < iter; ){  
    fips_inner_(S, Elist0, K, Sigma, Scc_list, Scc_inv_list);
    ++itcount;      
    
    if (true){ 
      switch(convcrit){
      case 1: 
	mad = mean_abs_diff_on_Emat_(S, Sigma, Emat0, 0);
	conv_check = mad;
	PRINT_CONV_CHECK1;	
	break;
      case 2:
	CONV_CHECK_LOGL_DIFF;
	break;
      case 3:
	CONV_CHECK_LOGL_DIFF_REF;
	PRINT_CONV_CHECK3;
	break;
      }	
      if (conv_check < eps) break;
    }      
  }
  
  logL = ips_logL_(S, K, nobs);  
  RETURN_VALUE;
}




  // double z;
  // int M = Haux.n_rows;
  // int N = Sigma.n_rows;
  // mat A = mat(N, N);
  // mat B = Sigma.cols(cc0);
  // for (size_t i=0; i<N; i++){
  //   for (size_t j=i; j<N; j++){
  //     z = 0;
  //     for (size_t k=0; k<M; k++){
  // 	for (size_t l=0; l<M; l++){
  // 	  z += B(i,k) * Haux(k,l) * B(j,l);
  // 	}
  //     }
  //     A(i,j) = z;
  //     A(j,i) = z;
  //   }
  // }
  // Sigma -= A;

  
  // mat V = Sigma.cols(cc0);
  // Sigma -= V * Haux * V.t();

// inline void fips_ggm_update_cc_parm_(const mat& S, const uvec& cc0, mat& K, mat& Sigma)
// {
  // mat Scc      = S.submat(cc0, cc0);
  // mat Scc_inv  = inv(Scc);
  
  // mat Kcc      = K.submat(cc0, cc0);	
  // mat Sigmacc  = Sigma.submat(cc0, cc0);
  
  // mat Kstar    = inv(Sigmacc);
  // mat Kupd     = Scc_inv + Kcc - Kstar;
  // mat Haux     = Kstar - Kstar * Scc * Kstar;  // Was Saux

  // K.submat(cc0, cc0) = Kupd;

  // Sigma -= Sigma.cols(cc0) * (Haux * Sigma.rows(cc0));

  // mat V = Sigma.cols(cc0);
  // Sigma -= V * Haux * V.t();
// }


// void fips_inner_(const mat& S, const List& Elist0,
		 // mat& K, mat& Sigma, List& Scc_list, List& Scc_inv_list)
// {
  // for (int i=0; i < Elist0.length(); ++i)
    // {
      // Rcout << "in loop" << endl;
      // uvec cc0 = Elist0[i];
      // fips_ggm_update_cc_parm_(S, cc0, K, Sigma);
    // }
// }
