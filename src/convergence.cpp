#include "RcppArmadillo.h"
#include "grips-utils.h"
#include "general-utils.h"
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

//[[Rcpp::export]] 
double califa_(const mat& S, const mat& Sigma, const mat& Sigmaold)
{
  double N = S.n_rows;
  double d = 0; //accu(abs(Sigma - Sigmaold)); // / (N * N);
  //double r = max(max(abs(Sigma)));
  double r = 0;
  for (int i=0; i<N; i++){
    for (int j=0; j<=i; j++){
      // double w = fabs(S(i,j));
      if (i != j)
	r += fabs(S(i,j));

      d += fabs(Sigma(i,j) - Sigmaold(i,j));
      // Rcout <<  " i : " << i << " j : " << j <<  " w ::: " << w << "  r ::: " << r << std::endl;
    }
  }
  r = r / (N * (N-1) / 2);
  d = d / (N * (N+1) / 2);
  double rat = d / r;

  // Rcout << "d : " << d << " r : " << r << " rat : " << rat << std::endl;
  
  return rat;
}



// --- END Utilities ---


// ----------------------------------------------------------------------
// Implements ||S(G) - Sigma(G)||_\infty
//
// Likelihood equation condition for GGM
// ----------------------------------------------------------------------


//[[Rcpp::export]] 
vec diff_on_Emat_(const mat& S, const mat& Sigma, const umat& E, int shift=1){
  uvec eids = sub2ind(size(S), E - shift);           // Obtain Element IDs
  vec v = S.elem(eids) - Sigma.elem(eids);   // Values of the Elements
  vec d = S.diag() - Sigma.diag();
  vec out = join_cols(d, v);
  return out;
}

//[[Rcpp::export]] 
vec diff_on_Elist_(const mat& S, const mat& Sigma, const List& E, int shift=1){
  umat Emat = conv_to<umat>::from(list2Emat_(E, shift));
  // Emat.print();
  return diff_on_Emat_(S, Sigma, Emat, 0);
  // return diff_on_Emat_(S, Sigma, Emat, shift);

  // uvec eids = sub2ind( size(S), Emat);           // Obtain Element IDs
  // vec v = S.elem( eids ) - Sigma.elem( eids );   // Values of the Elements
  // vec d = S.diag() - Sigma.diag();
  // vec out = join_cols(d, v);
  // return out;
  
}

//[[Rcpp::export]] 
double max_abs_diff_on_Emat_(const mat& Sigma, const mat& S, const umat& E, int shift=1){
  vec d = diff_on_Emat_(Sigma, S, E, shift);
  return max(abs(d));
}

//[[Rcpp::export]] 
double mean_abs_diff_on_Emat_(const mat& Sigma, const mat& S, const umat& E, int shift=1){
  vec d = diff_on_Emat_(Sigma, S, E, shift);
  return (sum(abs(d)) / d.size());
}

//[[Rcpp::export]] 
double max_abs_diff_on_Elist_(const mat& Sigma, const mat& S, const List& E, int shift=1){
  vec d = diff_on_Elist_(Sigma, S, E, shift);
  return max(abs(d));
}

//[[Rcpp::export]] 
double mean_abs_diff_on_Elist_(const mat& Sigma, const mat& S, const List& E, int shift=1){
  vec d = diff_on_Elist_(Sigma, S, E, shift);
  return (sum(abs(d)) / d.size());
}









// ----------------------------------------------------------------------
// implements ||{S(G) - Sigma(G)}^+||_\infty
//
// Likelihood equation condition for MTP2
// ----------------------------------------------------------------------

//[[Rcpp::export]] 
double max_diff_on_Emat_(const mat& Sigma, const mat& S, const mat& E){
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
// Implements || { S(G) - Sigma(G) } * K(G) ||_\infty
//
// Likelihood equation condition for MTP2
// ----------------------------------------------------------------------

//[[Rcpp::export]] 
double max_abs_diff_on_EK_(const mat& S, const mat& Sigma, const mat& E, const mat& K){
  double out=0, dd; 
  int u, v;
  
  for (size_t i=0; i<Sigma.n_rows; ++i)
    out = MAX(out, fabs(Sigma(i, i) - S(i, i)) * K(i, i));
 
  for (size_t j=0; j<E.n_cols; ++j){
    u = E(0,j) - 1; v = E(1,j) - 1 ;
    // v1 = Sigma(u, v);    v2 = S(u, v);    v3 = K(u, v);    dd  = fabs((v1 - v2) * v3);
    dd  = fabs((Sigma(u, v) - S(u, v)) * K(u, v));
    //Rcout << v1 << " " << v2 << " " << out << std::endl;
    out = MAX(out, dd);	    
  }
  return out;
}	






/* ********************************************************
 *
 * does edge fit to data
 *
 * ********************************************************/

// //[[Rcpp::export]]
bool does_edge_fit_to_data_ips(mat& S, uvec& uv, mat& K, mat& Sigma, double& eps){
  
  double d = 0;
  uvec uv0 = uv - 1;
  bool does_fit = true;  
  double maxvar = 1; //max_abs_diag_(S); // FIXME: Check only on uv

  //Rcout << "++ does_edge_fit_to_data_ips " << as_num(uv) << std::endl;

  mat dif = Sigma.submat(uv0, uv0) - S.submat(uv0, uv0);
  d = max_abs_(dif);
  
  does_fit = d < (eps * maxvar);
  return does_fit;  
}
  
// //[[Rcpp::export]]
bool does_edge_fit_to_data_mtp2(mat& S, uvec& uv, mat& K, mat& Sigma, double& eps){
  double d = 0, d01, dd;
  uvec uv0 = uv - 1;
  int u=uv0(0), v=uv0(1);
  bool does_fit=true;
  double maxvar = 1; //max_abs_diag_(S); // FIXME: Check only on uv

  // Rf_PrintValue(wrap(uv0));
  // Rcout << "u, v" << u << v << std::endl;
  // FIXME implement: d = max_abs_diag_diff_on_uv_(Sigma, S, uv0)
  mat dif = Sigma.submat(uv0, uv0) - S.submat(uv0, uv0);

  // Diagonal
  d = max_abs_diag_(dif);
  does_fit = d < (eps * maxvar);
  if (!(does_fit)) return false;

  // Off-diagonal  
  d01 = dif(0, 1);
  does_fit = d01 > -eps * maxvar;
  if (!(does_fit)) return false;
  
  dd = fabs(d01 * K(u, v));
  does_fit = dd < eps;
  if (!(does_fit)) return false;
  
  return does_fit;  
}	


// Rcout << " does edge fit mtp2= " << " does_fit=" << does_fit << std::endl;

// //[[Rcpp::export]]
bool does_edge_fit_to_data_lasso(mat& S, uvec& uv, mat& K, mat& Sigma, double& eps, mat& lambda){
  int u, v;
  bool does_fit=true;
  double duv, kuv, luv;
  
  for (size_t j=0; j<=1; ++j){
    for (size_t i=j; i<=1; ++i){
      u=uv(j) - 1;
      v=uv(i) - 1;
      duv = Sigma(u, v) - S(u, v);
      kuv = K(u, v);
      luv = lambda(u, v);

      if (u == v){
	if(erp>1)Rcout << "+++ does_edge_fit_to_data_lasso On diagonal" << std::endl;
	if (fabs(duv - luv) > eps) {
	  does_fit = false;
	}	
      } else { // On off-diagonal
	if (erp > 1)Rcout << "+++ does_edge_fit_to_data_lasso On off-diagonal" << std::endl;
	if (kuv >= eps){
	  // Rcout << "kuv >= eps" << std::endl;	
	  if (fabs(duv - luv) > eps) does_fit=false;
	} else {
	  if (kuv <= -eps){
	    //  Rcout << "kuv <= -eps" << std::endl;
	    if (fabs(duv + luv) > eps) does_fit=false;
	  }
	  else {
	    // Rcout << "kuv = 0" << std::endl;
	    if (fabs(duv) - luv > eps) does_fit=false;
	  }
	}
      }
      if(erp>=3) Rcout << "+++ does_edge_fit_to_data_lasso: " << as_num(uv) <<
		   " u = " << u << " v = " << v <<
		   " duv = " << duv << " kuv = " <<
		   kuv << " luv = " << luv << " does_fit = " << does_fit <<std::endl;
      if (!does_fit) goto exit;
    }
  }
  
 exit: {};
  return does_fit;
}

// //[[Rcpp::export]]
bool does_edge_fit_to_data_hybrid(mat& S, uvec& uv, mat& K, mat& Sigma, double& eps){
  Rcout << "does_edge_fit_to_data_hybrid (cpp) not implemented" << std::endl;
  return false;
}

// //[[Rcpp::export]]
bool does_edge_fit_to_data (mat& S, uvec& uv, mat& K, mat& Sigma, double& eps, mat& lambda, int imeth=0){
  
  switch(imeth){
  case 0: return does_edge_fit_to_data_ips   (S, uv, K, Sigma, eps); 
  case 1: return does_edge_fit_to_data_mtp2  (S, uv, K, Sigma, eps); 
  case 2: return does_edge_fit_to_data_lasso (S, uv, K, Sigma, eps, lambda); 
  case 3: return does_edge_fit_to_data_hybrid(S, uv, K, Sigma, eps);   
  default:
    Rcout << "does edge fit dunnowhattodo" << std::endl;    
  }
  return true;
}


/* ********************************************************
 *
 * does model fit to data
 *
 * ********************************************************/

// NOTE NOTE does_model_fit_to_data_ips computes average!!!
// //[[Rcpp::export]]
bool does_model_fit_to_data_ips(mat& S, umat& EE, mat& K, mat& Sigma,
				double& eps){
  // Rcout << "does_model_fit_to_data_ips (cpp)" << std::endl;
  double dd1, maxvar = 1; //max_abs_diag_(S);
  double N = S.n_rows;
  
  dd1 = max_abs_diff_on_Emat_(S, Sigma, EE) / N;
  bool does_fit = dd1 < (maxvar * eps);
  //printf(".does_model_fit_to_data_ips (cpp) d=%15e, s=%15e, does_fit=%i\n", d, s, does_fit);
  return does_fit;   
}

// //[[Rcpp::export]]
bool does_model_fit_to_data_mtp2(mat& S, umat& EE, mat& K, mat& Sigma,
				 double& eps){
  // Rcout << "does_model_fit_to_data_mtp2 (cpp)" << std::endl;
  int u, v;
  bool does_fit=true;
  double dd1, dd2, maxvar = 1; //max_abs_diag_(S);
  
  // Diagonal
  dd1 = max_abs_diag_diff_(S, Sigma);
  does_fit = dd1 < (maxvar * eps);
  if (!(does_fit)) return false;

  // Off-diagonal
  for (size_t j=0; j<EE.n_cols; j++){
    u = EE(0,j) - 1;
    v = EE(1,j) - 1 ;
    
    dd2 = Sigma(u, v) - S(u, v);
    does_fit = dd2 > -eps * maxvar;
    if (!(does_fit)) return false;

    dd1 = fabs(dd2 * K(u, v));
    does_fit = dd1 < eps;
    if (!(does_fit)) return false;
  }

  return does_fit;   
}


// //[[Rcpp::export]]
bool does_model_fit_to_data_lasso(mat& S, umat& EE, mat& K, mat& Sigma,
				  double& eps, mat& lambda){
  // Rcout << "does_model_fit_to_data_lasso (cpp)" << std::endl;
  int u, v;
  bool does_fit=true;
  double duv, kuv, luv;
  //  double maxvar = 1; //max_abs_diag_(S);

  //EE.print();
  for (size_t k=0; k<EE.n_cols; k++){
    vec uv = conv_to<vec>::from(EE.col(k));

    // Same code as for an edge ->
    for (size_t j=0; j<=1; ++j){
      for (size_t i=j; i<=1; ++i){
	// u = EE(0,j) - 1;
	// v = EE(1,j) - 1 ;
	u=uv(j) - 1;
	v=uv(i) - 1;
	duv = Sigma(u,v) - S(u,v);
	kuv = K(u,v);
	luv = lambda(u,v);
	
	if (u == v){
	  if(erp>1)Rcout << "+++ does_model_fit_to_data_lasso On diagonal" << std::endl;
	  if (fabs(duv - luv) > eps) {
	    does_fit=false;
	  }	
	} else { // On off-diagonal
	  if(erp>1)Rcout << "+++ does_model_fit_to_data_lasso On off-diagonal" << std::endl;
	  if (kuv >= eps){
	    // Rcout << "kuv >= eps" << std::endl;	
	    if (fabs(duv - luv) > eps) does_fit=false;
	  } else {
	    if (kuv <= -eps){
	      //  Rcout << "kuv <= -eps" << std::endl;
	      if (fabs(duv + luv) > eps) does_fit=false;
	    }
	    else {
	      // Rcout << "kuv = 0" << std::endl;
	      if (fabs(duv) - luv > eps) does_fit=false;
	    }
	  }
	}
	if(erp>=3) Rcout << "+++ does_model_fit_to_data_lasso: " << as_num(uv) <<
		     " u = " << u << " v = " << v <<
		     " duv = " << duv << " kuv = " <<
		     kuv << " luv = " << luv << " does_fit = " << does_fit <<std::endl;
	if (!does_fit) goto exit;
      }
    }
    
    // <- until here
    
  }
 exit: {};
  return does_fit;
}
  
// //[[Rcpp::export]]
bool does_model_fit_to_data_hybrid(mat& S, umat& EE, mat& K, mat& Sigma,
				   double& eps){
  Rcout << "does_model_fit_to_data_hybrid (cpp) not implemented" << std::endl;
  return false;   
}


// //[[Rcpp::export]]
bool does_model_fit_to_data(mat& S, umat& EE, mat& K, mat& Sigma,
			    double& eps, mat& lambda, int imeth=0){
    
  switch(imeth){
  case 0: return does_model_fit_to_data_ips   (S, EE, K, Sigma, eps);
  case 1: return does_model_fit_to_data_mtp2  (S, EE, K, Sigma, eps);	
  case 2: return does_model_fit_to_data_lasso (S, EE, K, Sigma, eps, lambda);	
  case 3: return does_model_fit_to_data_hybrid(S, EE, K, Sigma, eps);	
  default:
    Rcout << "does model fit dunnowhattodo" << std::endl;    
  }
  return true;
}


bool does_model_fit_to_data_relaxed(mat& S, mat& Sigma, mat& Sigma_prev, double& eps){

  double d2=0, d3=0;

  Rcout << "Sigma_prev (before)" << std::endl;  Sigma_prev.print();

  int dim = S.n_rows;
  
  //d2 = accu(abs(trimatl(S, -1))); // FIXME : constant across iterations

  d2 = accu(abs(S)); // FIXME : SLL is this ok???
  d2 = d2 / (dim * (dim - 1) / 2);

  //d3 = accu(abs(trimatl(Sigma) - trimatl(Sigma_prev))); 

  d3 = accu(abs(Sigma - Sigma_prev)); // FIXME : SLL is this ok???
  d3 = d3 / ((dim + 1) * dim / 2);
  
  Sigma_prev = Sigma;
  Rcout << "S " << std::endl;  S.print();  
  Rcout << "Sigma " << std::endl;  Sigma.print();  
  Rcout << "Sigma_prev (after)" << std::endl;  Sigma_prev.print();  
  Rcout << " d2 = " << d2 << " d2 * eps = " << d2 * eps <<" d3 = " << d3 << std::endl;

  bool does_fit = d3 < eps * d2;
  return does_fit;  
}

