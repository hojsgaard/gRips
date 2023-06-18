// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/gRips.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// inv_qr_
mat inv_qr_(mat& X);
RcppExport SEXP _gRips_inv_qr_(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(inv_qr_(X));
    return rcpp_result_gen;
END_RCPP
}
// setdiff_
arma::uvec setdiff_(arma::uvec x, arma::uvec y);
RcppExport SEXP _gRips_setdiff_(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(setdiff_(x, y));
    return rcpp_result_gen;
END_RCPP
}
// rep_nout
arma::vec rep_nout(vec x, unsigned int nout);
RcppExport SEXP _gRips_rep_nout(SEXP xSEXP, SEXP noutSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nout(noutSEXP);
    rcpp_result_gen = Rcpp::wrap(rep_nout(x, nout));
    return rcpp_result_gen;
END_RCPP
}
// vec2mat
arma::mat vec2mat(arma::vec x, int nrow, int ncol);
RcppExport SEXP _gRips_vec2mat(SEXP xSEXP, SEXP nrowSEXP, SEXP ncolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int >::type ncol(ncolSEXP);
    rcpp_result_gen = Rcpp::wrap(vec2mat(x, nrow, ncol));
    return rcpp_result_gen;
END_RCPP
}
// remove_elem
arma::vec remove_elem(arma::vec X, arma::uvec ent_, int shift);
RcppExport SEXP _gRips_remove_elem(SEXP XSEXP, SEXP ent_SEXP, SEXP shiftSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ent_(ent_SEXP);
    Rcpp::traits::input_parameter< int >::type shift(shiftSEXP);
    rcpp_result_gen = Rcpp::wrap(remove_elem(X, ent_, shift));
    return rcpp_result_gen;
END_RCPP
}
// extract_elem
arma::vec extract_elem(arma::vec X, arma::uvec ent_, int shift);
RcppExport SEXP _gRips_extract_elem(SEXP XSEXP, SEXP ent_SEXP, SEXP shiftSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ent_(ent_SEXP);
    Rcpp::traits::input_parameter< int >::type shift(shiftSEXP);
    rcpp_result_gen = Rcpp::wrap(extract_elem(X, ent_, shift));
    return rcpp_result_gen;
END_RCPP
}
// replace_elem
arma::vec replace_elem(arma::vec X, arma::uvec ent_, arma::vec value_, int shift);
RcppExport SEXP _gRips_replace_elem(SEXP XSEXP, SEXP ent_SEXP, SEXP value_SEXP, SEXP shiftSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ent_(ent_SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type value_(value_SEXP);
    Rcpp::traits::input_parameter< int >::type shift(shiftSEXP);
    rcpp_result_gen = Rcpp::wrap(replace_elem(X, ent_, value_, shift));
    return rcpp_result_gen;
END_RCPP
}
// remove_rows
arma::mat remove_rows(arma::mat X, arma::uvec row_ent_, int shift);
RcppExport SEXP _gRips_remove_rows(SEXP XSEXP, SEXP row_ent_SEXP, SEXP shiftSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type row_ent_(row_ent_SEXP);
    Rcpp::traits::input_parameter< int >::type shift(shiftSEXP);
    rcpp_result_gen = Rcpp::wrap(remove_rows(X, row_ent_, shift));
    return rcpp_result_gen;
END_RCPP
}
// extract_rows
arma::mat extract_rows(arma::mat X, arma::uvec row_ent_, int shift);
RcppExport SEXP _gRips_extract_rows(SEXP XSEXP, SEXP row_ent_SEXP, SEXP shiftSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type row_ent_(row_ent_SEXP);
    Rcpp::traits::input_parameter< int >::type shift(shiftSEXP);
    rcpp_result_gen = Rcpp::wrap(extract_rows(X, row_ent_, shift));
    return rcpp_result_gen;
END_RCPP
}
// replace_rows
arma::mat replace_rows(arma::mat X, arma::uvec row_ent_, arma::vec value, int shift);
RcppExport SEXP _gRips_replace_rows(SEXP XSEXP, SEXP row_ent_SEXP, SEXP valueSEXP, SEXP shiftSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type row_ent_(row_ent_SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type value(valueSEXP);
    Rcpp::traits::input_parameter< int >::type shift(shiftSEXP);
    rcpp_result_gen = Rcpp::wrap(replace_rows(X, row_ent_, value, shift));
    return rcpp_result_gen;
END_RCPP
}
// replace_u_vc
arma::mat replace_u_vc(arma::mat& M, arma::uvec u, arma::uvec v, arma::vec value, int shift);
RcppExport SEXP _gRips_replace_u_vc(SEXP MSEXP, SEXP uSEXP, SEXP vSEXP, SEXP valueSEXP, SEXP shiftSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type value(valueSEXP);
    Rcpp::traits::input_parameter< int >::type shift(shiftSEXP);
    rcpp_result_gen = Rcpp::wrap(replace_u_vc(M, u, v, value, shift));
    return rcpp_result_gen;
END_RCPP
}
// replace_uc_v
arma::mat replace_uc_v(arma::mat& M, arma::uvec u, arma::uvec v, arma::vec value, int shift);
RcppExport SEXP _gRips_replace_uc_v(SEXP MSEXP, SEXP uSEXP, SEXP vSEXP, SEXP valueSEXP, SEXP shiftSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type value(valueSEXP);
    Rcpp::traits::input_parameter< int >::type shift(shiftSEXP);
    rcpp_result_gen = Rcpp::wrap(replace_uc_v(M, u, v, value, shift));
    return rcpp_result_gen;
END_RCPP
}
// replace_uc_vc
arma::mat replace_uc_vc(arma::mat& M, arma::uvec u, arma::uvec v, arma::vec value, int shift);
RcppExport SEXP _gRips_replace_uc_vc(SEXP MSEXP, SEXP uSEXP, SEXP vSEXP, SEXP valueSEXP, SEXP shiftSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type value(valueSEXP);
    Rcpp::traits::input_parameter< int >::type shift(shiftSEXP);
    rcpp_result_gen = Rcpp::wrap(replace_uc_vc(M, u, v, value, shift));
    return rcpp_result_gen;
END_RCPP
}
// replace_u_v
arma::mat replace_u_v(arma::mat& M, arma::uvec u, arma::uvec v, arma::vec value, int shift);
RcppExport SEXP _gRips_replace_u_v(SEXP MSEXP, SEXP uSEXP, SEXP vSEXP, SEXP valueSEXP, SEXP shiftSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type value(valueSEXP);
    Rcpp::traits::input_parameter< int >::type shift(shiftSEXP);
    rcpp_result_gen = Rcpp::wrap(replace_u_v(M, u, v, value, shift));
    return rcpp_result_gen;
END_RCPP
}
// replace_uv_
arma::mat replace_uv_(arma::mat& M, arma::ivec u, arma::ivec v, arma::vec value, int shift);
RcppExport SEXP _gRips_replace_uv_(SEXP MSEXP, SEXP uSEXP, SEXP vSEXP, SEXP valueSEXP, SEXP shiftSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type value(valueSEXP);
    Rcpp::traits::input_parameter< int >::type shift(shiftSEXP);
    rcpp_result_gen = Rcpp::wrap(replace_uv_(M, u, v, value, shift));
    return rcpp_result_gen;
END_RCPP
}
// as_colvec
void as_colvec(arma::mat& M);
RcppExport SEXP _gRips_as_colvec(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type M(MSEXP);
    as_colvec(M);
    return R_NilValue;
END_RCPP
}
// extract_uv_
arma::mat extract_uv_(arma::mat& M, arma::ivec u, arma::ivec v, int shift);
RcppExport SEXP _gRips_extract_uv_(SEXP MSEXP, SEXP uSEXP, SEXP vSEXP, SEXP shiftSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type shift(shiftSEXP);
    rcpp_result_gen = Rcpp::wrap(extract_uv_(M, u, v, shift));
    return rcpp_result_gen;
END_RCPP
}
// make_complement_
int_vec make_complement_(int_vec cc, int d, int shift);
RcppExport SEXP _gRips_make_complement_(SEXP ccSEXP, SEXP dSEXP, SEXP shiftSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int_vec >::type cc(ccSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type shift(shiftSEXP);
    rcpp_result_gen = Rcpp::wrap(make_complement_(cc, d, shift));
    return rcpp_result_gen;
END_RCPP
}
// make_complement_list_
List make_complement_list_(List gen_lst, int d, int shift);
RcppExport SEXP _gRips_make_complement_list_(SEXP gen_lstSEXP, SEXP dSEXP, SEXP shiftSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type gen_lst(gen_lstSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type shift(shiftSEXP);
    rcpp_result_gen = Rcpp::wrap(make_complement_list_(gen_lst, d, shift));
    return rcpp_result_gen;
END_RCPP
}
// Scc_list_
List Scc_list_(const mat& S, const List& edges0);
RcppExport SEXP _gRips_Scc_list_(SEXP SSEXP, SEXP edges0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const List& >::type edges0(edges0SEXP);
    rcpp_result_gen = Rcpp::wrap(Scc_list_(S, edges0));
    return rcpp_result_gen;
END_RCPP
}
// Scc_inv_list_
List Scc_inv_list_(const mat& S, const List& edges0);
RcppExport SEXP _gRips_Scc_inv_list_(SEXP SSEXP, SEXP edges0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const List& >::type edges0(edges0SEXP);
    rcpp_result_gen = Rcpp::wrap(Scc_inv_list_(S, edges0));
    return rcpp_result_gen;
END_RCPP
}
// conips_ggm_
List conips_ggm_(arma::mat& S, List& elst, umat& emat, int& nobs, arma::mat K, int& maxit, double& eps, int& convcrit, int& print, List& aux);
RcppExport SEXP _gRips_conips_ggm_(SEXP SSEXP, SEXP elstSEXP, SEXP ematSEXP, SEXP nobsSEXP, SEXP KSEXP, SEXP maxitSEXP, SEXP epsSEXP, SEXP convcritSEXP, SEXP printSEXP, SEXP auxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< List& >::type elst(elstSEXP);
    Rcpp::traits::input_parameter< umat& >::type emat(ematSEXP);
    Rcpp::traits::input_parameter< int& >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K(KSEXP);
    Rcpp::traits::input_parameter< int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int& >::type convcrit(convcritSEXP);
    Rcpp::traits::input_parameter< int& >::type print(printSEXP);
    Rcpp::traits::input_parameter< List& >::type aux(auxSEXP);
    rcpp_result_gen = Rcpp::wrap(conips_ggm_(S, elst, emat, nobs, K, maxit, eps, convcrit, print, aux));
    return rcpp_result_gen;
END_RCPP
}
// covips_outer0_
List covips_outer0_(mat& S, mat& K, List& elst0, mat& Sigma, List& Scc_lst, List& Scci_lst, int& nobs, umat& emat_c, int& n_upd, double& max_visits, double& n_visits, double eps, int print);
RcppExport SEXP _gRips_covips_outer0_(SEXP SSEXP, SEXP KSEXP, SEXP elst0SEXP, SEXP SigmaSEXP, SEXP Scc_lstSEXP, SEXP Scci_lstSEXP, SEXP nobsSEXP, SEXP emat_cSEXP, SEXP n_updSEXP, SEXP max_visitsSEXP, SEXP n_visitsSEXP, SEXP epsSEXP, SEXP printSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< mat& >::type K(KSEXP);
    Rcpp::traits::input_parameter< List& >::type elst0(elst0SEXP);
    Rcpp::traits::input_parameter< mat& >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< List& >::type Scc_lst(Scc_lstSEXP);
    Rcpp::traits::input_parameter< List& >::type Scci_lst(Scci_lstSEXP);
    Rcpp::traits::input_parameter< int& >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< umat& >::type emat_c(emat_cSEXP);
    Rcpp::traits::input_parameter< int& >::type n_upd(n_updSEXP);
    Rcpp::traits::input_parameter< double& >::type max_visits(max_visitsSEXP);
    Rcpp::traits::input_parameter< double& >::type n_visits(n_visitsSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type print(printSEXP);
    rcpp_result_gen = Rcpp::wrap(covips_outer0_(S, K, elst0, Sigma, Scc_lst, Scci_lst, nobs, emat_c, n_upd, max_visits, n_visits, eps, print));
    return rcpp_result_gen;
END_RCPP
}
// covips_ggm_
List covips_ggm_(mat& S, List& elst, umat& emat, int& nobs, mat& K, int& maxit, double& eps, int& convcrit, int& print, List& aux);
RcppExport SEXP _gRips_covips_ggm_(SEXP SSEXP, SEXP elstSEXP, SEXP ematSEXP, SEXP nobsSEXP, SEXP KSEXP, SEXP maxitSEXP, SEXP epsSEXP, SEXP convcritSEXP, SEXP printSEXP, SEXP auxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< List& >::type elst(elstSEXP);
    Rcpp::traits::input_parameter< umat& >::type emat(ematSEXP);
    Rcpp::traits::input_parameter< int& >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< mat& >::type K(KSEXP);
    Rcpp::traits::input_parameter< int& >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int& >::type convcrit(convcritSEXP);
    Rcpp::traits::input_parameter< int& >::type print(printSEXP);
    Rcpp::traits::input_parameter< List& >::type aux(auxSEXP);
    rcpp_result_gen = Rcpp::wrap(covips_ggm_(S, elst, emat, nobs, K, maxit, eps, convcrit, print, aux));
    return rcpp_result_gen;
END_RCPP
}
// ncd_outer1_
List ncd_outer1_(mat& Sigma, mat& K, umat& emat, umat& emat_c, mat& amat, int& nobs, double& eps, int max_visits, int& n_visits, int print);
RcppExport SEXP _gRips_ncd_outer1_(SEXP SigmaSEXP, SEXP KSEXP, SEXP ematSEXP, SEXP emat_cSEXP, SEXP amatSEXP, SEXP nobsSEXP, SEXP epsSEXP, SEXP max_visitsSEXP, SEXP n_visitsSEXP, SEXP printSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< mat& >::type K(KSEXP);
    Rcpp::traits::input_parameter< umat& >::type emat(ematSEXP);
    Rcpp::traits::input_parameter< umat& >::type emat_c(emat_cSEXP);
    Rcpp::traits::input_parameter< mat& >::type amat(amatSEXP);
    Rcpp::traits::input_parameter< int& >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type max_visits(max_visitsSEXP);
    Rcpp::traits::input_parameter< int& >::type n_visits(n_visitsSEXP);
    Rcpp::traits::input_parameter< int >::type print(printSEXP);
    rcpp_result_gen = Rcpp::wrap(ncd_outer1_(Sigma, K, emat, emat_c, amat, nobs, eps, max_visits, n_visits, print));
    return rcpp_result_gen;
END_RCPP
}
// ncd_ggm_
List ncd_ggm_(mat& S, List& elst, umat& emat, int& nobs, mat K, int maxit, double& eps, int& convcrit, int print, List& aux);
RcppExport SEXP _gRips_ncd_ggm_(SEXP SSEXP, SEXP elstSEXP, SEXP ematSEXP, SEXP nobsSEXP, SEXP KSEXP, SEXP maxitSEXP, SEXP epsSEXP, SEXP convcritSEXP, SEXP printSEXP, SEXP auxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< List& >::type elst(elstSEXP);
    Rcpp::traits::input_parameter< umat& >::type emat(ematSEXP);
    Rcpp::traits::input_parameter< int& >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< mat >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int& >::type convcrit(convcritSEXP);
    Rcpp::traits::input_parameter< int >::type print(printSEXP);
    Rcpp::traits::input_parameter< List& >::type aux(auxSEXP);
    rcpp_result_gen = Rcpp::wrap(ncd_ggm_(S, elst, emat, nobs, K, maxit, eps, convcrit, print, aux));
    return rcpp_result_gen;
END_RCPP
}
// clone_
SEXP clone_(SEXP& x);
RcppExport SEXP _gRips_clone_(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(clone_(x));
    return rcpp_result_gen;
END_RCPP
}
// list_names_
chr_vec list_names_(List lst);
RcppExport SEXP _gRips_list_names_(SEXP lstSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type lst(lstSEXP);
    rcpp_result_gen = Rcpp::wrap(list_names_(lst));
    return rcpp_result_gen;
END_RCPP
}
// as_emat2amat_
mat as_emat2amat_(umat emat, int d);
RcppExport SEXP _gRips_as_emat2amat_(SEXP ematSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< umat >::type emat(ematSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(as_emat2amat_(emat, d));
    return rcpp_result_gen;
END_RCPP
}
// as_emat_complement_
umat as_emat_complement_(umat emat, int d);
RcppExport SEXP _gRips_as_emat_complement_(SEXP ematSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< umat >::type emat(ematSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(as_emat_complement_(emat, d));
    return rcpp_result_gen;
END_RCPP
}
// has_full_rank_
bool has_full_rank_(mat& Delta, double eps);
RcppExport SEXP _gRips_has_full_rank_(SEXP DeltaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type Delta(DeltaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(has_full_rank_(Delta, eps));
    return rcpp_result_gen;
END_RCPP
}
// project_onto_G_
mat project_onto_G_(const mat& Delta, const umat& emc);
RcppExport SEXP _gRips_project_onto_G_(SEXP DeltaSEXP, SEXP emcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type Delta(DeltaSEXP);
    Rcpp::traits::input_parameter< const umat& >::type emc(emcSEXP);
    rcpp_result_gen = Rcpp::wrap(project_onto_G_(Delta, emc));
    return rcpp_result_gen;
END_RCPP
}
// mnorm_one_
double mnorm_one_(mat& Delta);
RcppExport SEXP _gRips_mnorm_one_(SEXP DeltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type Delta(DeltaSEXP);
    rcpp_result_gen = Rcpp::wrap(mnorm_one_(Delta));
    return rcpp_result_gen;
END_RCPP
}
// mnorm_maxabs_
double mnorm_maxabs_(mat& Delta);
RcppExport SEXP _gRips_mnorm_maxabs_(SEXP DeltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type Delta(DeltaSEXP);
    rcpp_result_gen = Rcpp::wrap(mnorm_maxabs_(Delta));
    return rcpp_result_gen;
END_RCPP
}
// ggm_logL_
double ggm_logL_(mat& S, mat& K, int nobs);
RcppExport SEXP _gRips_ggm_logL_(SEXP SSEXP, SEXP KSEXP, SEXP nobsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< mat& >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    rcpp_result_gen = Rcpp::wrap(ggm_logL_(S, K, nobs));
    return rcpp_result_gen;
END_RCPP
}

// validate (ensure exported C++ functions exist before calling them)
static int _gRips_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _gRips_RcppExport_registerCCallable() { 
    R_RegisterCCallable("gRips", "_gRips_RcppExport_validate", (DL_FUNC)_gRips_RcppExport_validate);
    return R_NilValue;
}
