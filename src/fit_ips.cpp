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
List conips_ggm_(arma::mat& S, List& elst, umat& emat, int& nobs,
	      arma::mat K,
	      int& maxit, double& eps, int& convcrit, int& print, List& aux){

  int version      = aux["version"];  
  double logL, logLp, conv_check=9999, gap=-1.0;
//  double conv_ref; // FIXME needed??
  int itcount  = 0;
  double nparm = S.n_cols + emat.n_cols;  
  umat emat0   = emat - 1;
  umat emat_c  = as_emat_complement_(emat-1, S.n_rows);
  double mno, maxabs;
  mat Sigma, dif, Delta;

  List elst0 = clone(elst);
  for (int i=0; i<elst.length(); i++){
    elst0[i] = as<arma::uvec>(elst0[i]) - 1; // 0-based
  }
  
  // Variables for CONIPS only
  List clist0 = make_clist_(S, elst);  
  for (int i=0; i<elst.length(); i++){
    clist0[i] = as<arma::uvec>(clist0[i]) - 1; // 0-based			     
  }  
  // END - Variables for CONIPS only  

  
  INIT_CONVERGENCE_CHECK;
    
  for (; itcount < maxit; ){  
    conips_inner_(S, K, elst0, clist0, print=print);
    ++itcount;          

    switch(convcrit){
    case 1: 
      Sigma = inv_qr_(K); // Needed for conips only
      dif   = S - Sigma;
      Delta = project_onto_G_(dif, emat_c);
      maxabs = mnorm_maxabs_(Delta);
      conv_check = maxabs;
      if (print>=3){
	mno   = mnorm_one_(Delta);
	logL  = ggm_logL_(S, K, nobs);
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

  Sigma = inv_qr_(K); // FOR CONIPS ONLY

  logL  = ggm_logL_(S, K, nobs);
  RETURN_VALUE;
}

  

// ------------------------------------------------------------------
// ------ COVIPS - Fast iterative proportional scaling          -------
// ------------------------------------------------------------------



void covips_update_parm0_(const uvec& cc0, const mat& Scc, mat& K, mat& Sigma,
			  const mat& Scc_inv,
			  int& nupd, double eps=0.01, int print=0)
{
  
  mat Sigmacc  = Sigma.submat(cc0, cc0);
  mat Kcc, Kupd, Haux, B, dd;

  dd = Scc - Sigmacc;
  double Knorm = mnorm_maxabs_(dd);
  
  if (Knorm > eps){
    if (print>=4)
      Rprintf(">>>> K norm %f - yes do update \n", Knorm);
    mat Kstar = inv(Sigmacc);
    Kcc      = K.submat(cc0, cc0);	
    Kupd     = Scc_inv + Kcc - Kstar;
    Haux     = Kstar - Kstar * Scc * Kstar;  

    K.submat(cc0, cc0) = Kupd;
    mat B  = Sigma.cols(cc0);
    Sigma -= B * Haux * B.t();
    nupd++;
  } else {
    if (print>=4)    
      Rprintf(">>>> K norm %f - no do not update \n", Knorm);
  }
}


void covips_update_parm_(const uvec& cc0, const mat& Scc, mat& K, mat& Sigma,
			 const mat& Scc_inv, int print=0)
{  
  mat Sigmacc  = Sigma.submat(cc0, cc0);
  mat Kcc, Kupd, Haux;
  mat Kstar    = inv(Sigmacc);
  
  Kcc      = K.submat(cc0, cc0);	
  Kupd     = Scc_inv + Kcc - Kstar;
  Haux     = Kstar - Kstar * Scc * Kstar;  

  K.submat(cc0, cc0)  = Kupd;
  mat B = Sigma.cols(cc0);
  Sigma -= B * Haux * B.t();  
}


// Update some edges once
void covips_inner0_(mat& S, mat& K, List& elst0,
		    mat& Sigma, List& Scc_lst, List& Scci_lst,
		    int& nupd, double eps=0.01, int print=0)
{
  nupd = 0;
  for (int i=0; i < elst0.length(); ++i){
    uvec cc0 = elst0[i];
    covips_update_parm0_(cc0, Scc_lst[i], K, Sigma, Scci_lst[i],
			 nupd=nupd, eps=eps, print=print);
  }
}

// Update all edges once
List covips_inner_(const mat& S, mat& K, const List& elst0,
		   mat& Sigma, List& Scc_lst, List& Scci_lst,
		   int print=0)
{
  for (int i=0; i < elst0.length(); ++i){
    uvec cc0 = elst0[i];
    covips_update_parm_(cc0, Scc_lst[i], K, Sigma, Scci_lst[i],
			       print=print);
  }
  return List::create(_["iter"] = 1);
}



//[[Rcpp::export]]  
List covips_outer0_(mat& S, mat& K, List& elst0, mat& Sigma,
		   List& Scc_lst, List& Scci_lst,
		   int& nobs, int& maxit,
		   int& nupd, double eps=0.01, int print=0){

  int itcount = 0;
  double logL;
  
  for (; itcount < maxit; ){  
    covips_inner0_(S=S, K=K, elst0=elst0, Sigma=Sigma,
  		  Scc_lst=Scc_lst, Scci_lst=Scci_lst,
  		  nupd=nupd, eps=eps, print=print);
    if (print>=3){
      logL = ggm_logL_(S, K, nobs);  
      Rprintf(">>> itcount: %3d maxit: %d nupd: %6d edges: %4d logL: %16.8f\n",
	      itcount, maxit, nupd, elst0.length(), logL);
    }
    ++itcount;
    if (nupd <= 0){
      break;        
    }  
  }
  return List::create(_["iter"] = itcount);
}

//[[Rcpp::export(.c_covips_ggm_)]] 
List covips_ggm_(mat& S, List& elst, umat& emat, int& nobs,
	       mat& K,       
	       int& maxit, double& eps, int& convcrit, int& print, List& aux){

  int version      = aux["version"];
  double logL, logLp, conv_check=9999, gap=-1.0, mno, maxabs;
  double nparm = S.n_cols + emat.n_cols;
  umat emat0   = emat - 1;  
  umat emat_c  = as_emat_complement_(emat-1, S.n_rows);
  mat Sigma, dif, Delta;  
  List res1, res2;
  int iter1=0, iter2=0, itcount=0, nupd=0;
  
  List elst0 = clone(elst);  
  for (int i=0; i < elst.length(); i++) {
    elst0[i] = as<arma::uvec>(elst0[i]) - 1; // 0-based
  }
    
  // Variables for COVIPS only
  Sigma = inv_qr_(K);
    
  List Scc_lst  = Scc_list_(S, elst0);  
  List Scci_lst = Scc_inv_list_(S, elst0);
  // END - Variables for COVIPS only  

  INIT_CONVERGENCE_CHECK;
  
  res1 = covips_outer0_(S=S, K=K, elst0=elst0, Sigma=Sigma,
			Scc_lst=Scc_lst, Scci_lst=Scci_lst,
			nobs=nobs, maxit=maxit, nupd=nupd,
			eps=eps, print=print);
  iter1 = res1["iter"];

  dif    = S - Sigma;
  Delta  = project_onto_G_(dif, emat_c);
  maxabs = mnorm_maxabs_(Delta);
  conv_check = maxabs;

  if (print>=2)
    Rprintf(">> iterations for start: %d maxabs: %14.10f\n", iter1, maxabs);

   version=1;
  if (version==1){
    while (true){
      res2 = covips_inner_(S=S, K=K, elst0=elst0, Sigma=Sigma,
			   Scc_lst=Scc_lst, Scci_lst=Scci_lst,
			   print=print);
      ++iter2;
      
      switch(convcrit){
      case 1: 	
	dif    = S - Sigma;
	Delta  = project_onto_G_(dif, emat_c);
	maxabs = mnorm_maxabs_(Delta);
	conv_check = maxabs;
	if (print >= 3){
	  mno  = mnorm_one_(Delta);
	  logL = ggm_logL_(S, K, nobs);  
	  Rprintf(">>> covips iter: %4d eps: %14.10f, mno: %14.10f maxabs: %14.10f logL: %14.10f\n",
		  iter2, eps, mno, maxabs, logL);
	}
	break;
      case 2:
	CONV_CHECK_LOGL_DIFF;
	break;
      }
      if ((iter2 == maxit) || (conv_check < eps)) break;
    }
  }
  
  itcount = iter2 + iter1;
  logL    = ggm_logL_(S, K, nobs);  
  RETURN_VALUE;
}

