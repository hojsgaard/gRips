#include "RcppArmadillo.h"
using namespace Rcpp;
using namespace arma;

typedef Rcpp::NumericVector   num_vec;
typedef Rcpp::IntegerVector   int_vec;
typedef Rcpp::CharacterVector chr_vec;

const int erp=0; // ERROR PRINTING

// const double MEPS = 1e-6;

#define as_num(cc) NumericVector(cc.begin(), cc.end())

mat project_K_onto_G_(const mat& K, const umat& emc);
double mnormone_(mat& Delta);
int method2int_(CharacterVector method);
double ggm_logL_(mat& S, mat& K, int nobs);
mat initSigma_(mat& S);
mat initK_(mat& S);
double get_conv_ref(const List& aux);



// -------------------------------------------------------------
// Macros common to COVIPS and CONIPS
// -------------------------------------------------------------

#define RETURN_VALUE					\
  return List::create(					\
  _["K"]     = K,					\
  _["Sigma"] = Sigma,					\
  _["logL"]  = logL,					\
  _["iter"]  = itcount,					\
  _["gap"]   = gap,					\
  _["conv_check"] = conv_check);			\


#define INIT_CONVERGENCE_CHECK			\
  switch (convcrit){                            \
  case 1:					\
    break;					\
  case 2:					\
    logLp = ggm_logL_(S, K, nobs);		\
    break;					\
  case 3:					\
    conv_ref = get_conv_ref(aux);		\
    break;                                      \
  case 4:                                       \
    break;                                      \
  }                                             \

#define PRINT_CONV_CHECK1			                \
  if ((print > 0) && ((itcount % 1) == 0)) {			\
  logL  = ggm_logL_(S, K, nobs);				\
  sprintf(buffer, "itcount=%4i, logL=%f, nparm=%4f, mad=%f\n",	\
	  itcount, logL, nparm, mad);				\
  Rcout << buffer;						\
  }								\


#define PRINT_CONV_CHECK3						\
  if ((print > 0) && ((itcount % 1) == 0)) {				\
    sprintf(buffer, "itcount=%4i, logL=%f, logLp=%f, nparm=%4f, conv_ref=%f, conv_check=%f\n", \
	    itcount, logL, logLp, nparm, conv_ref, conv_check);		\
    Rcout << buffer;							\
  }									\

#define CONV_CHECK_LOGL_DIFF			\
  logL  = ggm_logL_(S, K, nobs);		\
  conv_check = fabs(logL - logLp)  / nparm;	\
  logLp = logL;					\
  
#define CONV_CHECK_LOGL_DIFF_REF			\
  logL  = ggm_logL_(S, K, nobs);			\
  conv_check = fabs(logL - conv_ref) / nparm;		\
  logLp = logL;						\

