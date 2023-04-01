#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _gRips_RcppExport_registerCCallable();
extern SEXP _gRips_Scc_inv_list_(SEXP, SEXP);
extern SEXP _gRips_Scc_list_(SEXP, SEXP);
extern SEXP _gRips_Sigma_to_K_(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_as_colvec(SEXP);
extern SEXP _gRips_as_emat2amat_(SEXP, SEXP);
extern SEXP _gRips_as_emat_complement_(SEXP, SEXP);
extern SEXP _gRips_califa_(SEXP, SEXP, SEXP);
extern SEXP _gRips_clone_(SEXP);
extern SEXP _gRips_conips_ggm_(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_covips_ggm_(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_coxips_ggm_(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_diff_fun_(SEXP, SEXP, SEXP);
extern SEXP _gRips_diff_on_Elist_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_diff_on_emat_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_extract_elem(SEXP, SEXP, SEXP);
extern SEXP _gRips_extract_rows(SEXP, SEXP, SEXP);
extern SEXP _gRips_extract_uv_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_find_str_(SEXP, SEXP);
extern SEXP _gRips_ggm_logL_(SEXP, SEXP, SEXP);
extern SEXP _gRips_inv_qr_(SEXP);
extern SEXP _gRips_list2col_(SEXP, SEXP);
extern SEXP _gRips_list2row_(SEXP, SEXP);
extern SEXP _gRips_list_names_(SEXP);
extern SEXP _gRips_list_to_emat(SEXP, SEXP);
extern SEXP _gRips_make_clist_(SEXP, SEXP);
extern SEXP _gRips_max_abs_(SEXP);
extern SEXP _gRips_max_abs_diag_(SEXP);
extern SEXP _gRips_max_abs_diag_diff_(SEXP, SEXP);
extern SEXP _gRips_max_abs_diff_(SEXP, SEXP);
extern SEXP _gRips_max_abs_diff_on_EK_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_max_abs_diff_on_Elist_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_max_abs_diff_on_emat_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_max_abs_diff_rel_(SEXP, SEXP);
extern SEXP _gRips_max_diag_diff_(SEXP, SEXP);
extern SEXP _gRips_max_diff_on_emat_(SEXP, SEXP, SEXP);
extern SEXP _gRips_mean_abs_diff_on_Elist_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_mean_abs_diff_on_emat_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_mnormone_(SEXP);
extern SEXP _gRips_ncd_ggm_(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_outerloop1_(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_project_onto_G_(SEXP, SEXP);
extern SEXP _gRips_remove_elem(SEXP, SEXP, SEXP);
extern SEXP _gRips_remove_rows(SEXP, SEXP, SEXP);
extern SEXP _gRips_rep_nout(SEXP, SEXP);
extern SEXP _gRips_replace_elem(SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_replace_rows(SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_replace_u_v(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_replace_u_vc(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_replace_uc_v(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_replace_uc_vc(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_replace_uv_(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRips_setdiff_(SEXP, SEXP);
extern SEXP _gRips_unique_cols(SEXP);
extern SEXP _gRips_unique_rows(SEXP);
extern SEXP _gRips_vec2mat(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_gRips_RcppExport_registerCCallable", (DL_FUNC) &_gRips_RcppExport_registerCCallable,  0},
    {"_gRips_Scc_inv_list_",                (DL_FUNC) &_gRips_Scc_inv_list_,                 2},
    {"_gRips_Scc_list_",                    (DL_FUNC) &_gRips_Scc_list_,                     2},
    {"_gRips_Sigma_to_K_",                  (DL_FUNC) &_gRips_Sigma_to_K_,                   5},
    {"_gRips_as_colvec",                    (DL_FUNC) &_gRips_as_colvec,                     1},
    {"_gRips_as_emat2amat_",                (DL_FUNC) &_gRips_as_emat2amat_,                 2},
    {"_gRips_as_emat_complement_",          (DL_FUNC) &_gRips_as_emat_complement_,           2},
    {"_gRips_califa_",                      (DL_FUNC) &_gRips_califa_,                       3},
    {"_gRips_clone_",                       (DL_FUNC) &_gRips_clone_,                        1},
    {"_gRips_conips_ggm_",                  (DL_FUNC) &_gRips_conips_ggm_,                  10},
    {"_gRips_covips_ggm_",                  (DL_FUNC) &_gRips_covips_ggm_,                  10},
    {"_gRips_coxips_ggm_",                  (DL_FUNC) &_gRips_coxips_ggm_,                  10},
    {"_gRips_diff_fun_",                    (DL_FUNC) &_gRips_diff_fun_,                     3},
    {"_gRips_diff_on_Elist_",               (DL_FUNC) &_gRips_diff_on_Elist_,                4},
    {"_gRips_diff_on_emat_",                (DL_FUNC) &_gRips_diff_on_emat_,                 4},
    {"_gRips_extract_elem",                 (DL_FUNC) &_gRips_extract_elem,                  3},
    {"_gRips_extract_rows",                 (DL_FUNC) &_gRips_extract_rows,                  3},
    {"_gRips_extract_uv_",                  (DL_FUNC) &_gRips_extract_uv_,                   4},
    {"_gRips_find_str_",                    (DL_FUNC) &_gRips_find_str_,                     2},
    {"_gRips_ggm_logL_",                    (DL_FUNC) &_gRips_ggm_logL_,                     3},
    {"_gRips_inv_qr_",                      (DL_FUNC) &_gRips_inv_qr_,                       1},
    {"_gRips_list2col_",                    (DL_FUNC) &_gRips_list2col_,                     2},
    {"_gRips_list2row_",                    (DL_FUNC) &_gRips_list2row_,                     2},
    {"_gRips_list_names_",                  (DL_FUNC) &_gRips_list_names_,                   1},
    {"_gRips_list_to_emat",                 (DL_FUNC) &_gRips_list_to_emat,                  2},
    {"_gRips_make_clist_",                  (DL_FUNC) &_gRips_make_clist_,                   2},
    {"_gRips_max_abs_",                     (DL_FUNC) &_gRips_max_abs_,                      1},
    {"_gRips_max_abs_diag_",                (DL_FUNC) &_gRips_max_abs_diag_,                 1},
    {"_gRips_max_abs_diag_diff_",           (DL_FUNC) &_gRips_max_abs_diag_diff_,            2},
    {"_gRips_max_abs_diff_",                (DL_FUNC) &_gRips_max_abs_diff_,                 2},
    {"_gRips_max_abs_diff_on_EK_",          (DL_FUNC) &_gRips_max_abs_diff_on_EK_,           4},
    {"_gRips_max_abs_diff_on_Elist_",       (DL_FUNC) &_gRips_max_abs_diff_on_Elist_,        4},
    {"_gRips_max_abs_diff_on_emat_",        (DL_FUNC) &_gRips_max_abs_diff_on_emat_,         4},
    {"_gRips_max_abs_diff_rel_",            (DL_FUNC) &_gRips_max_abs_diff_rel_,             2},
    {"_gRips_max_diag_diff_",               (DL_FUNC) &_gRips_max_diag_diff_,                2},
    {"_gRips_max_diff_on_emat_",            (DL_FUNC) &_gRips_max_diff_on_emat_,             3},
    {"_gRips_mean_abs_diff_on_Elist_",      (DL_FUNC) &_gRips_mean_abs_diff_on_Elist_,       4},
    {"_gRips_mean_abs_diff_on_emat_",       (DL_FUNC) &_gRips_mean_abs_diff_on_emat_,        4},
    {"_gRips_mnormone_",                    (DL_FUNC) &_gRips_mnormone_,                     1},
    {"_gRips_ncd_ggm_",                     (DL_FUNC) &_gRips_ncd_ggm_,                     10},
    {"_gRips_outerloop1_",                  (DL_FUNC) &_gRips_outerloop1_,                   9},
    {"_gRips_project_onto_G_",            (DL_FUNC) &_gRips_project_onto_G_,             2},
    {"_gRips_remove_elem",                  (DL_FUNC) &_gRips_remove_elem,                   3},
    {"_gRips_remove_rows",                  (DL_FUNC) &_gRips_remove_rows,                   3},
    {"_gRips_rep_nout",                     (DL_FUNC) &_gRips_rep_nout,                      2},
    {"_gRips_replace_elem",                 (DL_FUNC) &_gRips_replace_elem,                  4},
    {"_gRips_replace_rows",                 (DL_FUNC) &_gRips_replace_rows,                  4},
    {"_gRips_replace_u_v",                  (DL_FUNC) &_gRips_replace_u_v,                   5},
    {"_gRips_replace_u_vc",                 (DL_FUNC) &_gRips_replace_u_vc,                  5},
    {"_gRips_replace_uc_v",                 (DL_FUNC) &_gRips_replace_uc_v,                  5},
    {"_gRips_replace_uc_vc",                (DL_FUNC) &_gRips_replace_uc_vc,                 5},
    {"_gRips_replace_uv_",                  (DL_FUNC) &_gRips_replace_uv_,                   5},
    {"_gRips_setdiff_",                     (DL_FUNC) &_gRips_setdiff_,                      2},
    {"_gRips_unique_cols",                  (DL_FUNC) &_gRips_unique_cols,                   1},
    {"_gRips_unique_rows",                  (DL_FUNC) &_gRips_unique_rows,                   1},
    {"_gRips_vec2mat",                      (DL_FUNC) &_gRips_vec2mat,                       3},
    {NULL, NULL, 0}
};

void R_init_gRips(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
