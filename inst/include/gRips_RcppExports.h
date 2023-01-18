// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_gRips_RCPPEXPORTS_H_GEN_
#define RCPP_gRips_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace gRips {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("gRips", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("gRips", "_gRips_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in gRips");
            }
        }
    }

    inline double max_abs_(const mat& S) {
        typedef SEXP(*Ptr_max_abs_)(SEXP);
        static Ptr_max_abs_ p_max_abs_ = NULL;
        if (p_max_abs_ == NULL) {
            validateSignature("double(*max_abs_)(const mat&)");
            p_max_abs_ = (Ptr_max_abs_)R_GetCCallable("gRips", "_gRips_max_abs_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_max_abs_(Shield<SEXP>(Rcpp::wrap(S)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double max_abs_diag_(const mat& S) {
        typedef SEXP(*Ptr_max_abs_diag_)(SEXP);
        static Ptr_max_abs_diag_ p_max_abs_diag_ = NULL;
        if (p_max_abs_diag_ == NULL) {
            validateSignature("double(*max_abs_diag_)(const mat&)");
            p_max_abs_diag_ = (Ptr_max_abs_diag_)R_GetCCallable("gRips", "_gRips_max_abs_diag_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_max_abs_diag_(Shield<SEXP>(Rcpp::wrap(S)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double max_abs_diff_rel_(const mat& S, const mat& Sigma) {
        typedef SEXP(*Ptr_max_abs_diff_rel_)(SEXP,SEXP);
        static Ptr_max_abs_diff_rel_ p_max_abs_diff_rel_ = NULL;
        if (p_max_abs_diff_rel_ == NULL) {
            validateSignature("double(*max_abs_diff_rel_)(const mat&,const mat&)");
            p_max_abs_diff_rel_ = (Ptr_max_abs_diff_rel_)R_GetCCallable("gRips", "_gRips_max_abs_diff_rel_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_max_abs_diff_rel_(Shield<SEXP>(Rcpp::wrap(S)), Shield<SEXP>(Rcpp::wrap(Sigma)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double max_abs_diff_(const mat& S, const mat& Sigma) {
        typedef SEXP(*Ptr_max_abs_diff_)(SEXP,SEXP);
        static Ptr_max_abs_diff_ p_max_abs_diff_ = NULL;
        if (p_max_abs_diff_ == NULL) {
            validateSignature("double(*max_abs_diff_)(const mat&,const mat&)");
            p_max_abs_diff_ = (Ptr_max_abs_diff_)R_GetCCallable("gRips", "_gRips_max_abs_diff_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_max_abs_diff_(Shield<SEXP>(Rcpp::wrap(S)), Shield<SEXP>(Rcpp::wrap(Sigma)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double max_abs_diag_diff_(const mat& S, const mat& Sigma) {
        typedef SEXP(*Ptr_max_abs_diag_diff_)(SEXP,SEXP);
        static Ptr_max_abs_diag_diff_ p_max_abs_diag_diff_ = NULL;
        if (p_max_abs_diag_diff_ == NULL) {
            validateSignature("double(*max_abs_diag_diff_)(const mat&,const mat&)");
            p_max_abs_diag_diff_ = (Ptr_max_abs_diag_diff_)R_GetCCallable("gRips", "_gRips_max_abs_diag_diff_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_max_abs_diag_diff_(Shield<SEXP>(Rcpp::wrap(S)), Shield<SEXP>(Rcpp::wrap(Sigma)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double max_diag_diff_(const mat& S, const mat& Sigma) {
        typedef SEXP(*Ptr_max_diag_diff_)(SEXP,SEXP);
        static Ptr_max_diag_diff_ p_max_diag_diff_ = NULL;
        if (p_max_diag_diff_ == NULL) {
            validateSignature("double(*max_diag_diff_)(const mat&,const mat&)");
            p_max_diag_diff_ = (Ptr_max_diag_diff_)R_GetCCallable("gRips", "_gRips_max_diag_diff_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_max_diag_diff_(Shield<SEXP>(Rcpp::wrap(S)), Shield<SEXP>(Rcpp::wrap(Sigma)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double califa_(const mat& S, const mat& Sigma, const mat& Sigmaold) {
        typedef SEXP(*Ptr_califa_)(SEXP,SEXP,SEXP);
        static Ptr_califa_ p_califa_ = NULL;
        if (p_califa_ == NULL) {
            validateSignature("double(*califa_)(const mat&,const mat&,const mat&)");
            p_califa_ = (Ptr_califa_)R_GetCCallable("gRips", "_gRips_califa_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_califa_(Shield<SEXP>(Rcpp::wrap(S)), Shield<SEXP>(Rcpp::wrap(Sigma)), Shield<SEXP>(Rcpp::wrap(Sigmaold)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline vec diff_on_Emat_(const mat& S, const mat& Sigma, const umat& E, int shift = 1) {
        typedef SEXP(*Ptr_diff_on_Emat_)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_diff_on_Emat_ p_diff_on_Emat_ = NULL;
        if (p_diff_on_Emat_ == NULL) {
            validateSignature("vec(*diff_on_Emat_)(const mat&,const mat&,const umat&,int)");
            p_diff_on_Emat_ = (Ptr_diff_on_Emat_)R_GetCCallable("gRips", "_gRips_diff_on_Emat_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_diff_on_Emat_(Shield<SEXP>(Rcpp::wrap(S)), Shield<SEXP>(Rcpp::wrap(Sigma)), Shield<SEXP>(Rcpp::wrap(E)), Shield<SEXP>(Rcpp::wrap(shift)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<vec >(rcpp_result_gen);
    }

    inline vec diff_on_Elist_(const mat& S, const mat& Sigma, const List& E, int shift = 1) {
        typedef SEXP(*Ptr_diff_on_Elist_)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_diff_on_Elist_ p_diff_on_Elist_ = NULL;
        if (p_diff_on_Elist_ == NULL) {
            validateSignature("vec(*diff_on_Elist_)(const mat&,const mat&,const List&,int)");
            p_diff_on_Elist_ = (Ptr_diff_on_Elist_)R_GetCCallable("gRips", "_gRips_diff_on_Elist_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_diff_on_Elist_(Shield<SEXP>(Rcpp::wrap(S)), Shield<SEXP>(Rcpp::wrap(Sigma)), Shield<SEXP>(Rcpp::wrap(E)), Shield<SEXP>(Rcpp::wrap(shift)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<vec >(rcpp_result_gen);
    }

    inline double max_abs_diff_on_Emat_(const mat& Sigma, const mat& S, const umat& E, int shift = 1) {
        typedef SEXP(*Ptr_max_abs_diff_on_Emat_)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_max_abs_diff_on_Emat_ p_max_abs_diff_on_Emat_ = NULL;
        if (p_max_abs_diff_on_Emat_ == NULL) {
            validateSignature("double(*max_abs_diff_on_Emat_)(const mat&,const mat&,const umat&,int)");
            p_max_abs_diff_on_Emat_ = (Ptr_max_abs_diff_on_Emat_)R_GetCCallable("gRips", "_gRips_max_abs_diff_on_Emat_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_max_abs_diff_on_Emat_(Shield<SEXP>(Rcpp::wrap(Sigma)), Shield<SEXP>(Rcpp::wrap(S)), Shield<SEXP>(Rcpp::wrap(E)), Shield<SEXP>(Rcpp::wrap(shift)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double mean_abs_diff_on_Emat_(const mat& Sigma, const mat& S, const umat& E, int shift = 1) {
        typedef SEXP(*Ptr_mean_abs_diff_on_Emat_)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_mean_abs_diff_on_Emat_ p_mean_abs_diff_on_Emat_ = NULL;
        if (p_mean_abs_diff_on_Emat_ == NULL) {
            validateSignature("double(*mean_abs_diff_on_Emat_)(const mat&,const mat&,const umat&,int)");
            p_mean_abs_diff_on_Emat_ = (Ptr_mean_abs_diff_on_Emat_)R_GetCCallable("gRips", "_gRips_mean_abs_diff_on_Emat_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_mean_abs_diff_on_Emat_(Shield<SEXP>(Rcpp::wrap(Sigma)), Shield<SEXP>(Rcpp::wrap(S)), Shield<SEXP>(Rcpp::wrap(E)), Shield<SEXP>(Rcpp::wrap(shift)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double max_abs_diff_on_Elist_(const mat& Sigma, const mat& S, const List& E, int shift = 1) {
        typedef SEXP(*Ptr_max_abs_diff_on_Elist_)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_max_abs_diff_on_Elist_ p_max_abs_diff_on_Elist_ = NULL;
        if (p_max_abs_diff_on_Elist_ == NULL) {
            validateSignature("double(*max_abs_diff_on_Elist_)(const mat&,const mat&,const List&,int)");
            p_max_abs_diff_on_Elist_ = (Ptr_max_abs_diff_on_Elist_)R_GetCCallable("gRips", "_gRips_max_abs_diff_on_Elist_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_max_abs_diff_on_Elist_(Shield<SEXP>(Rcpp::wrap(Sigma)), Shield<SEXP>(Rcpp::wrap(S)), Shield<SEXP>(Rcpp::wrap(E)), Shield<SEXP>(Rcpp::wrap(shift)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double mean_abs_diff_on_Elist_(const mat& Sigma, const mat& S, const List& E, int shift = 1) {
        typedef SEXP(*Ptr_mean_abs_diff_on_Elist_)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_mean_abs_diff_on_Elist_ p_mean_abs_diff_on_Elist_ = NULL;
        if (p_mean_abs_diff_on_Elist_ == NULL) {
            validateSignature("double(*mean_abs_diff_on_Elist_)(const mat&,const mat&,const List&,int)");
            p_mean_abs_diff_on_Elist_ = (Ptr_mean_abs_diff_on_Elist_)R_GetCCallable("gRips", "_gRips_mean_abs_diff_on_Elist_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_mean_abs_diff_on_Elist_(Shield<SEXP>(Rcpp::wrap(Sigma)), Shield<SEXP>(Rcpp::wrap(S)), Shield<SEXP>(Rcpp::wrap(E)), Shield<SEXP>(Rcpp::wrap(shift)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double max_diff_on_Emat_(const mat& Sigma, const mat& S, const mat& E) {
        typedef SEXP(*Ptr_max_diff_on_Emat_)(SEXP,SEXP,SEXP);
        static Ptr_max_diff_on_Emat_ p_max_diff_on_Emat_ = NULL;
        if (p_max_diff_on_Emat_ == NULL) {
            validateSignature("double(*max_diff_on_Emat_)(const mat&,const mat&,const mat&)");
            p_max_diff_on_Emat_ = (Ptr_max_diff_on_Emat_)R_GetCCallable("gRips", "_gRips_max_diff_on_Emat_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_max_diff_on_Emat_(Shield<SEXP>(Rcpp::wrap(Sigma)), Shield<SEXP>(Rcpp::wrap(S)), Shield<SEXP>(Rcpp::wrap(E)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double max_abs_diff_on_EK_(const mat& S, const mat& Sigma, const mat& E, const mat& K) {
        typedef SEXP(*Ptr_max_abs_diff_on_EK_)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_max_abs_diff_on_EK_ p_max_abs_diff_on_EK_ = NULL;
        if (p_max_abs_diff_on_EK_ == NULL) {
            validateSignature("double(*max_abs_diff_on_EK_)(const mat&,const mat&,const mat&,const mat&)");
            p_max_abs_diff_on_EK_ = (Ptr_max_abs_diff_on_EK_)R_GetCCallable("gRips", "_gRips_max_abs_diff_on_EK_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_max_abs_diff_on_EK_(Shield<SEXP>(Rcpp::wrap(S)), Shield<SEXP>(Rcpp::wrap(Sigma)), Shield<SEXP>(Rcpp::wrap(E)), Shield<SEXP>(Rcpp::wrap(K)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline Rcpp::NumericMatrix list2row_(Rcpp::List input_list, int d = 2) {
        typedef SEXP(*Ptr_list2row_)(SEXP,SEXP);
        static Ptr_list2row_ p_list2row_ = NULL;
        if (p_list2row_ == NULL) {
            validateSignature("Rcpp::NumericMatrix(*list2row_)(Rcpp::List,int)");
            p_list2row_ = (Ptr_list2row_)R_GetCCallable("gRips", "_gRips_list2row_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_list2row_(Shield<SEXP>(Rcpp::wrap(input_list)), Shield<SEXP>(Rcpp::wrap(d)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::NumericMatrix >(rcpp_result_gen);
    }

    inline Rcpp::NumericMatrix list2col_(Rcpp::List input_list, int d = 2) {
        typedef SEXP(*Ptr_list2col_)(SEXP,SEXP);
        static Ptr_list2col_ p_list2col_ = NULL;
        if (p_list2col_ == NULL) {
            validateSignature("Rcpp::NumericMatrix(*list2col_)(Rcpp::List,int)");
            p_list2col_ = (Ptr_list2col_)R_GetCCallable("gRips", "_gRips_list2col_");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_list2col_(Shield<SEXP>(Rcpp::wrap(input_list)), Shield<SEXP>(Rcpp::wrap(d)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::NumericMatrix >(rcpp_result_gen);
    }

}

#endif // RCPP_gRips_RCPPEXPORTS_H_GEN_
