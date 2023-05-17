#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _gRips_RcppExport_registerCCallable(void);
extern SEXP _gRips_Scc_inv_list_(void *, void *);
extern SEXP _gRips_Scc_list_(void *, void *);
extern SEXP _gRips_as_colvec(void *);
extern SEXP _gRips_as_emat2amat_(void *, void *);
extern SEXP _gRips_as_emat_complement_(void *, void *);
extern SEXP _gRips_clone_(void *);
extern SEXP _gRips_conips_ggm_(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _gRips_covips_ggm_(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _gRips_covips_loop0_(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _gRips_diff_on_elst_(void *, void *, void *, void *);
extern SEXP _gRips_diff_on_emat_(void *, void *, void *, void *);
extern SEXP _gRips_extract_elem(void *, void *, void *);
extern SEXP _gRips_extract_rows(void *, void *, void *);
extern SEXP _gRips_extract_uv_(void *, void *, void *, void *);
extern SEXP _gRips_find_str_(void *, void *);
extern SEXP _gRips_ggm_logL_(void *, void *, void *);
extern SEXP _gRips_has_full_rank_(void *);
extern SEXP _gRips_inv_qr_(void *);
extern SEXP _gRips_list2col_(void *, void *);
extern SEXP _gRips_list2row_(void *, void *);
extern SEXP _gRips_list_names_(void *);
extern SEXP _gRips_list_to_emat(void *, void *);
extern SEXP _gRips_make_clist_(void *, void *);
extern SEXP _gRips_max_abs_(void *);
extern SEXP _gRips_max_abs_diag_(void *);
extern SEXP _gRips_max_abs_diag_diff_(void *, void *);
extern SEXP _gRips_max_abs_diff_(void *, void *);
extern SEXP _gRips_max_abs_diff_on_elst_(void *, void *, void *, void *);
extern SEXP _gRips_max_abs_diff_on_emat_(void *, void *, void *, void *);
extern SEXP _gRips_max_abs_diff_rel_(void *, void *);
extern SEXP _gRips_max_diag_diff_(void *, void *);
extern SEXP _gRips_max_diff_on_emat_(void *, void *, void *);
extern SEXP _gRips_mean_abs_diff_on_elst_(void *, void *, void *, void *);
extern SEXP _gRips_mean_abs_diff_on_emat_(void *, void *, void *, void *);
extern SEXP _gRips_mnorm_maxabs_(void *);
extern SEXP _gRips_mnorm_one_(void *);
extern SEXP _gRips_ncd_ggm_(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _gRips_outerloop1_(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP _gRips_project_onto_G_(void *, void *);
extern SEXP _gRips_remove_elem(void *, void *, void *);
extern SEXP _gRips_remove_rows(void *, void *, void *);
extern SEXP _gRips_rep_nout(void *, void *);
extern SEXP _gRips_replace_elem(void *, void *, void *, void *);
extern SEXP _gRips_replace_rows(void *, void *, void *, void *);
extern SEXP _gRips_replace_u_v(void *, void *, void *, void *, void *);
extern SEXP _gRips_replace_u_vc(void *, void *, void *, void *, void *);
extern SEXP _gRips_replace_uc_v(void *, void *, void *, void *, void *);
extern SEXP _gRips_replace_uc_vc(void *, void *, void *, void *, void *);
extern SEXP _gRips_replace_uv_(void *, void *, void *, void *, void *);
extern SEXP _gRips_setdiff_(void *, void *);
extern SEXP _gRips_unique_cols(void *);
extern SEXP _gRips_unique_rows(void *);
extern SEXP _gRips_vec2mat(void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"_gRips_RcppExport_registerCCallable", (DL_FUNC) &_gRips_RcppExport_registerCCallable,  0},
    {"_gRips_Scc_inv_list_",                (DL_FUNC) &_gRips_Scc_inv_list_,                 2},
    {"_gRips_Scc_list_",                    (DL_FUNC) &_gRips_Scc_list_,                     2},
    {"_gRips_as_colvec",                    (DL_FUNC) &_gRips_as_colvec,                     1},
    {"_gRips_as_emat2amat_",                (DL_FUNC) &_gRips_as_emat2amat_,                 2},
    {"_gRips_as_emat_complement_",          (DL_FUNC) &_gRips_as_emat_complement_,           2},
    {"_gRips_clone_",                       (DL_FUNC) &_gRips_clone_,                        1},
    {"_gRips_conips_ggm_",                  (DL_FUNC) &_gRips_conips_ggm_,                  10},
    {"_gRips_covips_ggm_",                  (DL_FUNC) &_gRips_covips_ggm_,                  10},
    {"_gRips_covips_loop0_",                (DL_FUNC) &_gRips_covips_loop0_,                11},
    {"_gRips_diff_on_elst_",                (DL_FUNC) &_gRips_diff_on_elst_,                 4},
    {"_gRips_diff_on_emat_",                (DL_FUNC) &_gRips_diff_on_emat_,                 4},
    {"_gRips_extract_elem",                 (DL_FUNC) &_gRips_extract_elem,                  3},
    {"_gRips_extract_rows",                 (DL_FUNC) &_gRips_extract_rows,                  3},
    {"_gRips_extract_uv_",                  (DL_FUNC) &_gRips_extract_uv_,                   4},
    {"_gRips_find_str_",                    (DL_FUNC) &_gRips_find_str_,                     2},
    {"_gRips_ggm_logL_",                    (DL_FUNC) &_gRips_ggm_logL_,                     3},
    {"_gRips_has_full_rank_",               (DL_FUNC) &_gRips_has_full_rank_,                1},
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
    {"_gRips_max_abs_diff_on_elst_",        (DL_FUNC) &_gRips_max_abs_diff_on_elst_,         4},
    {"_gRips_max_abs_diff_on_emat_",        (DL_FUNC) &_gRips_max_abs_diff_on_emat_,         4},
    {"_gRips_max_abs_diff_rel_",            (DL_FUNC) &_gRips_max_abs_diff_rel_,             2},
    {"_gRips_max_diag_diff_",               (DL_FUNC) &_gRips_max_diag_diff_,                2},
    {"_gRips_max_diff_on_emat_",            (DL_FUNC) &_gRips_max_diff_on_emat_,             3},
    {"_gRips_mean_abs_diff_on_elst_",       (DL_FUNC) &_gRips_mean_abs_diff_on_elst_,        4},
    {"_gRips_mean_abs_diff_on_emat_",       (DL_FUNC) &_gRips_mean_abs_diff_on_emat_,        4},
    {"_gRips_mnorm_maxabs_",                (DL_FUNC) &_gRips_mnorm_maxabs_,                 1},
    {"_gRips_mnorm_one_",                   (DL_FUNC) &_gRips_mnorm_one_,                    1},
    {"_gRips_ncd_ggm_",                     (DL_FUNC) &_gRips_ncd_ggm_,                     10},
    {"_gRips_outerloop1_",                  (DL_FUNC) &_gRips_outerloop1_,                   9},
    {"_gRips_project_onto_G_",              (DL_FUNC) &_gRips_project_onto_G_,               2},
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
