#include "RcppArmadillo.h"
#include "grips_utils.h"
#include "general_utils.h"
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

//[[Rcpp::depends(RcppClock)]]

#include <RcppClock.h>
#include <thread>

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)>(b))?(b):(a))


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
List Scc_list_(const mat& S, const List& edges0){
  List out(edges0.length());
  for (int i=0; i<edges0.length(); i++){
    uvec cc = as<arma::uvec>(edges0[i]);
    mat Scc = S.submat(cc, cc); 
    out[i] = Scc;
  }
  return out;
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


inline void conips_ggm_update_cc_parm_(mat& S, mat& K, uvec& cc0, uvec& aa0, mat& Scc_inv, int print=0)
{
  // mat Scc=S.submat(cc0, cc0);
  mat Kaa=K.submat(aa0, aa0);
  mat Kca=K.submat(cc0, aa0);    
  mat Kac=K.submat(aa0, cc0);    
  // mat L = Kca * inv(Kaa) * Kac;
  K.submat(cc0, cc0) = Scc_inv + Kca * inv_qr_(Kaa) * Kac;
}

void conips_inner_(arma::mat& S, arma::mat& K,
		   List& edges0, List& clist0, int print=0){  // NOTE: edges0, clist0 are 0-based

  List Scc_inv_list = Scc_inv_list_(S, edges0);

  // Rcout << "conips_inner started" << std::endl;
  for (int i=0; i < edges0.length(); i++)
    {    
      uvec cc0 = (edges0[i]), aa0 = (clist0[i]);
      mat Scc_inv = Scc_inv_list[i];
      conips_ggm_update_cc_parm_(S, K, cc0, aa0, Scc_inv, print);
    } 
  // Rcout << "conips_inner done" << std::endl;
}


//[[Rcpp::export(.c_conips_ggm_)]]
List conips_ggm_(arma::mat& S, List& elist, umat& emat, int& nobs,
	      arma::mat K,
	      int& maxiter, double& eps, int& convcrit, int& print, List& aux){

  double logL, logLp, conv_check=9999, gap=-1.0;
  double conv_ref; // FIXME needed??
  // char buffer[200];
  int itcount  = 0;
  double nparm = S.n_cols + emat.n_cols;  
  umat emat0   = emat - 1;

  umat emat_c  = as_emat_complement_(emat-1, S.n_rows);
  // double eps2  = MIN(eps, 1.0/S.n_rows);  
  mat dif, Delta;
  double mno, maxabs;

  // Variables for conips
  mat Sigma;
  List elist0 = clone(elist);
  List clist0 = make_clist_(S, elist);
  
  for (int i=0; i<elist.length(); i++){
    elist0[i] = as<arma::uvec>(elist0[i]) - 1; // 0-based
    clist0[i] = as<arma::uvec>(clist0[i]) - 1; // 0-based			     
  }
  // END
  
  INIT_CONVERGENCE_CHECK;
    
  for (; itcount < maxiter; ){  
    conips_inner_(S, K, elist0, clist0, print=print);
    ++itcount;          

    switch(convcrit){
    case 1: 
      Sigma = inv_qr_(K); // Needed for conips only
      // mad   = mean_abs_diff_on_emat_(S, Sigma, emat0, 0);
      // conv_check = mad;
      dif = S - Sigma;
      Delta = project_onto_G_(dif, emat_c);
      mno = mnorm_one_(Delta);
      conv_check = mno;
      if (print>=3)
	Rprintf(">>> conips iter: %4d eps: %14.10f, mno: %14.10f\n", itcount, eps, mno);
      // PRINT_CONV_CHECK1;
      break;
    case 2:
      CONV_CHECK_LOGL_DIFF;
      break;
    }    
    if (conv_check < eps) break;
  }

  Sigma = inv_qr_(K); // FOR CONIPS ONLY

  logL  = ggm_logL_(S, K, nobs);
  RETURN_VALUE;
}

  

// ------------------------------------------------------------------
// ------ COVIPS - Fast iterative proportional scaling          -------
// ------------------------------------------------------------------

void covips_ggm_update_cc_parm_(const mat& Scc, const uvec& cc0, mat& K, mat& Sigma,
				const mat& Scc_inv, int print=0)
{  
  mat Kcc      = K.submat(cc0, cc0);	
  mat Sigmacc  = Sigma.submat(cc0, cc0);
  mat Kstar    = inv(Sigmacc);
  mat Kupd     = Scc_inv + Kcc - Kstar;
  mat Haux     = Kstar - Kstar * Scc * Kstar;  // Was Saux
  K.submat(cc0, cc0)  = Kupd;
  Sigma              -= Sigma.cols(cc0) * Haux * Sigma.rows(cc0);
}

List covips_inner_(const mat& S, mat& K, 
		   const List& elist0,
		   mat& Sigma, List& Scc_lst, List& Scci_lst,
		   int& nupdates,
		   int smart=0, double eps_smart=0.0, int print=0)
{
  for (int i=0; i < elist0.length(); ++i){
    uvec cc0 = elist0[i];
    covips_ggm_update_cc_parm_(Scc_lst[i], cc0, K, Sigma, Scci_lst[i],
			       print=print);
  }
  return List::create(_["iter"] = 1);
}

void covips_ggm_update_cc_parm0_(const mat& Scc, const uvec& cc0, mat& K, mat& Sigma,
				 const mat& Scc_inv,
				 int& nupdates,
				 int smart=0, double eps_smart=0.0, int print=0)
{
  
  mat Sigmacc  = Sigma.submat(cc0, cc0);
  mat Kstar    = inv(Sigmacc);
  mat dd = Scc_inv - Kstar;
  double Knorm = mnorm_one_(dd);

  mat Kcc, Kupd, Haux;
  
  if (Knorm > 0.01){
    Kcc      = K.submat(cc0, cc0);	
    Kupd     = Scc_inv + Kcc - Kstar;
    Haux     = Kstar - Kstar * Scc * Kstar;  // Was Saux
    // cc0.t().print();
    if (print>=4)
      Rprintf(">>>> K norm %f - yes do update \n", Knorm);
    K.submat(cc0, cc0) = Kupd;
    Sigma -= Sigma.cols(cc0) * Haux * Sigma.rows(cc0);
    nupdates++;
  } else {
    // cc0.t().print();
    if (print>=4)    
      Rprintf(">>>> K norm %f - no do not update \n", Knorm);
  }
}

void covips_inner0_(mat& S, mat& K, 
		    List& elist0,
		    mat& Sigma, List& Scc_lst, List& Scci_lst,
		    int& nupdates,
		    int smart=0, double eps_smart=0.0, int print=0)
{
  nupdates = 0;
  for (int i=0; i < elist0.length(); ++i){
    uvec cc0 = elist0[i];
    covips_ggm_update_cc_parm0_(Scc_lst[i], cc0, K, Sigma, Scci_lst[i],
				nupdates=nupdates,
				print=print);
  }
}

//[[Rcpp::export]]  
List covips_loop0_(mat& S, mat& K, List& elist0, mat& Sigma,
		   List& Scc_lst, List& Scci_lst,
		   int& maxiter,
		   int nupdates, int smart=0, double eps_smart=0.0, int print=0){
  int itcount = 0;

  for (; itcount < maxiter; ){  
    covips_inner0_(S=S, K=K, elist0=elist0, Sigma=Sigma,
  		  Scc_lst=Scc_lst, Scci_lst=Scci_lst,
  		  nupdates=nupdates, smart=smart, eps_smart=eps_smart,
  		  print=print);
    if (print>=3)
      Rprintf(">>> itcount: %d maxiter: %d nupdates: %d\n", itcount, maxiter, nupdates);
    ++itcount;
    if (nupdates <= elist0.length() / 10)
      break;  
  }
  // Rprintf("itcount: %d\n", itcount);
  return List::create(_["iter"] = itcount);
}


// https://gallery.rcpp.org/articles/RcppClock-benchmarking-Rcpp-code/
  // Rcpp::Clock clock;
  
  // clock.tick("both_naps");
  // clock.tock("both_naps");
  
  // // send the times to the R global environment variable, named "naptimes"
  // clock.stop("naptimes");

//[[Rcpp::export(.c_covips_ggm_)]] 
List covips_ggm_(mat& S, List& elist, umat& emat, int& nobs,
	       mat& K,       
	       int& maxiter, double& eps, int& convcrit, int& print, List& aux){

  int version      = aux["version"];
  int  smart       = aux["smart"];
  double eps_smart = aux["eps_smart"];
  
  double logL, logLp, conv_check=9999, gap=-1.0;
  double conv_ref; // FIXME needed??
  int itcount  = 0;
  double nparm = S.n_cols + emat.n_cols;
  umat emat0   = emat - 1;  
  umat emat_c  = as_emat_complement_(emat-1, S.n_rows);
  mat dif, Delta;
  double mno, maxabs;
  mat Sigma    = initSigma_(S);  
  List res1, res2;
  int iter1=0, iter2=0;

  
  INIT_CONVERGENCE_CHECK;

  // Compute Scc etc only once
  List elist0 = clone(elist);  
  for (int i=0; i < elist.length(); i++) {
    elist0[i] = as<arma::uvec>(elist0[i]) - 1; // 0-based
  }
  
  List Scci_lst = Scc_inv_list_(S, elist0);
  List Scc_lst  = Scc_list_(S, elist0);
  int nupdates=0;
  
  if (version==0){
    res1 = covips_loop0_(S=S, K=K, elist0=elist0, Sigma=Sigma,
			 Scc_lst=Scc_lst, Scci_lst=Scci_lst,
			 maxiter=maxiter, nupdates=nupdates,
			 smart=smart, eps_smart=eps_smart, print=print);
    iter1 = res1["iter"];
    if (print>=2)
      Rprintf(">> iterations for start: %d \n", iter1);
  }

  
  for (; itcount < maxiter; ){  
    res2 = covips_inner_(S=S, K=K, elist0=elist0, Sigma=Sigma,
			 Scc_lst=Scc_lst, Scci_lst=Scci_lst,
			 nupdates=nupdates, smart=smart, eps_smart=eps_smart,
			 print=print);
    ++itcount;      
    switch(convcrit){
    case 1: 	
      // mad = mean_abs_diff_on_emat_(S, Sigma, emat0, 0);
      // conv_check = mad;
      dif = S - Sigma;
      Delta = project_onto_G_(dif, emat_c);
      mno = mnorm_one_(Delta);
      conv_check = mno;
      if (print>=3){
	maxabs = mnorm_maxabs_(Delta);
	logL = ggm_logL_(S, K, nobs);  
	Rprintf(">>> covips iter: %4d eps: %14.10f, mno: %14.10f maxabs: %14.10f logL: %14.10f\n", itcount, eps, mno, maxabs, logL);
      }
	// PRINT_CONV_CHECK1;	
      break;
    case 2:
      CONV_CHECK_LOGL_DIFF;
      break;
    }	
    if (conv_check < eps) break;
  }
  
  logL = ggm_logL_(S, K, nobs);  
  RETURN_VALUE;
}





//[[Rcpp::export(.c_coxips_ggm_)]] 
List coxips_ggm_(mat& S, List& elist, umat& emat, int& nobs,
		 mat K,       
		 int& maxiter, double& eps, int& convcrit, int& print, List& aux){

  Rprintf("coxips_ggm_\n");
  
  int version      = aux["version"];
  int smart        = aux["smart"];
  double eps_smart = aux["eps_smart"];

  mat Sigma;
  List elist0, clist0, Scci_lst, Scc_lst;

  double logL, logLp, conv_check=9999, gap=-1.0;
  // char buffer[200];
  int itcount  = 0;
  double nparm = S.n_cols + emat.n_cols;
  umat emat0   = emat - 1;
  
  umat emat_c  = as_emat_complement_(emat-1, S.n_rows);
  mat dif, Delta;
  double mno;
  int nupdates=0;
  
  // Initialization
  switch(version){
  case 1: // covips
    Sigma    = initSigma_(S);  
    elist0 = clone(elist);  
    for (int i=0; i < elist.length(); i++) {
      elist0[i] = as<arma::uvec>(elist0[i]) - 1; // 0-based
    }
    
    Scci_lst = Scc_inv_list_(S, elist0);
    Scc_lst  = Scc_list_(S, elist0);
    break;
  case 2: // conips
    elist0 = clone(elist);
    clist0 = make_clist_(S, elist);
  
    for (int i=0; i<elist.length(); i++){
      elist0[i] = as<arma::uvec>(elist0[i]) - 1; // 0-based
      clist0[i] = as<arma::uvec>(clist0[i]) - 1; // 0-based			     
    }

    break;
  }

  // Iteration
  for (; itcount < maxiter; ){  
    switch(version){
    case 1: // covips
      covips_inner_(S=S, K=K, elist0=elist0, Sigma=Sigma,
		    Scc_lst=Scc_lst, Scci_lst=Scci_lst,
		    nupdates=nupdates, smart=smart, eps_smart=eps_smart,
		    print=print);      
      break;
    case 2: // conips
      conips_inner_(S, K, elist0, clist0, print=print);
      break;
    }
    ++itcount;          

    switch(convcrit){
    case 1: 	
      Sigma = inv_qr_(K); // Needed for conips only
      dif = S - Sigma;
      Delta = project_onto_G_(dif, emat_c);
      mno = mnorm_one_(Delta);
      conv_check = mno;
      if (print>=3)
	Rprintf(">>> covips iter: %4d eps: %14.10f, mno: %14.10f\n", itcount, eps, mno);
      break;
    case 2:
      CONV_CHECK_LOGL_DIFF;
      break;
    }	
      
    if (conv_check < eps) break;
  }

  if (version == 2) // conips
    Sigma = inv_qr_(K); // FOR CONIPS ONLY

  logL = ggm_logL_(S, K, nobs);  
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




// void covips_ggm_update_cc_parm0_(const mat& Scc, const uvec& cc0, mat& K, mat& Sigma,
// 				      const mat& Scc_inv)
// {
  
//   mat Kcc      = K.submat(cc0, cc0);	
//   mat Sigmacc  = Sigma.submat(cc0, cc0);
  
//   mat Kstar  = inv(Sigmacc);
//   mat Kupd   = Scc_inv + Kcc - Kstar;
//   mat Haux   = Kstar - Kstar * Scc * Kstar;  // Was Saux

//   K.submat(cc0, cc0) = Kupd;

//   Sigma -= Sigma.cols(cc0) * Haux * Sigma.rows(cc0);
//   // mat V = Sigma.cols(cc0);
//   // Sigma -= V * Haux * V.t();
// }


  // if (smart == 1){
  //   nupdates++;
  //   Sigma -= Sigma.cols(cc0) * Haux * Sigma.rows(cc0);
  // } else {

  //   mat AAt = Sigma.rows(cc0) * Sigma.cols(cc0);
  //   mat HAAt = Haux * AAt ;
  //   double dd = accu(HAAt * HAAt) / accu(AAt); 
    
  //   if (dd > eps_smart){
  //     if (print>=4){
  // 	Rprintf(">>>> updating dd %12.9f", dd);cc0.t().print();
  //     }
  //     nupdates++;
  //     Sigma -= Sigma.cols(cc0) * Haux * Sigma.rows(cc0);     
  //   }
  //   else
  //     {
  // 	if (print>=4) {
  // 	  Rprintf("not updating dd %12.9f", dd);cc0.t().print();
  // 	}
  //     }
  // }
  // Sigma -= Sigma.cols(cc0) * (Kstar - Kstar * Scc * Kstar) * Sigma.rows(cc0);

