// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// logf_z_given_A
double logf_z_given_A(NumericVector zs, List A, NumericVector pis, NumericMatrix Ps);
RcppExport SEXP _mlsbm_logf_z_given_A(SEXP zsSEXP, SEXP ASEXP, SEXP pisSEXP, SEXP PsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type zs(zsSEXP);
    Rcpp::traits::input_parameter< List >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pis(pisSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Ps(PsSEXP);
    rcpp_result_gen = Rcpp::wrap(logf_z_given_A(zs, A, pis, Ps));
    return rcpp_result_gen;
END_RCPP
}
// logf_z_given_A_single
double logf_z_given_A_single(NumericVector zs, NumericMatrix A, NumericVector pis, NumericMatrix Ps);
RcppExport SEXP _mlsbm_logf_z_given_A_single(SEXP zsSEXP, SEXP ASEXP, SEXP pisSEXP, SEXP PsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type zs(zsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pis(pisSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Ps(PsSEXP);
    rcpp_result_gen = Rcpp::wrap(logf_z_given_A_single(zs, A, pis, Ps));
    return rcpp_result_gen;
END_RCPP
}
// rdirichlet_cpp
arma::mat rdirichlet_cpp(int num_samples, arma::vec alpha_m);
RcppExport SEXP _mlsbm_rdirichlet_cpp(SEXP num_samplesSEXP, SEXP alpha_mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type num_samples(num_samplesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha_m(alpha_mSEXP);
    rcpp_result_gen = Rcpp::wrap(rdirichlet_cpp(num_samples, alpha_m));
    return rcpp_result_gen;
END_RCPP
}
// sample_SBM_fast
NumericVector sample_SBM_fast(NumericVector z, NumericMatrix P);
RcppExport SEXP _mlsbm_sample_SBM_fast(SEXP zSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_SBM_fast(z, P));
    return rcpp_result_gen;
END_RCPP
}
// update_P
NumericMatrix update_P(List A, NumericVector zs, int K0, double b10, double b20);
RcppExport SEXP _mlsbm_update_P(SEXP ASEXP, SEXP zsSEXP, SEXP K0SEXP, SEXP b10SEXP, SEXP b20SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type zs(zsSEXP);
    Rcpp::traits::input_parameter< int >::type K0(K0SEXP);
    Rcpp::traits::input_parameter< double >::type b10(b10SEXP);
    Rcpp::traits::input_parameter< double >::type b20(b20SEXP);
    rcpp_result_gen = Rcpp::wrap(update_P(A, zs, K0, b10, b20));
    return rcpp_result_gen;
END_RCPP
}
// update_P_single
NumericMatrix update_P_single(NumericMatrix A, NumericVector zs, int K0, double b10, double b20);
RcppExport SEXP _mlsbm_update_P_single(SEXP ASEXP, SEXP zsSEXP, SEXP K0SEXP, SEXP b10SEXP, SEXP b20SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type zs(zsSEXP);
    Rcpp::traits::input_parameter< int >::type K0(K0SEXP);
    Rcpp::traits::input_parameter< double >::type b10(b10SEXP);
    Rcpp::traits::input_parameter< double >::type b20(b20SEXP);
    rcpp_result_gen = Rcpp::wrap(update_P_single(A, zs, K0, b10, b20));
    return rcpp_result_gen;
END_RCPP
}
// update_counts
NumericVector update_counts(NumericVector zs, int K0);
RcppExport SEXP _mlsbm_update_counts(SEXP zsSEXP, SEXP K0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type zs(zsSEXP);
    Rcpp::traits::input_parameter< int >::type K0(K0SEXP);
    rcpp_result_gen = Rcpp::wrap(update_counts(zs, K0));
    return rcpp_result_gen;
END_RCPP
}
// update_z
NumericVector update_z(NumericVector zs, List A, NumericMatrix Ps, NumericVector pis, NumericVector classes);
RcppExport SEXP _mlsbm_update_z(SEXP zsSEXP, SEXP ASEXP, SEXP PsSEXP, SEXP pisSEXP, SEXP classesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type zs(zsSEXP);
    Rcpp::traits::input_parameter< List >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Ps(PsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pis(pisSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type classes(classesSEXP);
    rcpp_result_gen = Rcpp::wrap(update_z(zs, A, Ps, pis, classes));
    return rcpp_result_gen;
END_RCPP
}
// update_z_single
NumericVector update_z_single(NumericVector zs, NumericMatrix A, NumericMatrix Ps, NumericVector pis, NumericVector classes);
RcppExport SEXP _mlsbm_update_z_single(SEXP zsSEXP, SEXP ASEXP, SEXP PsSEXP, SEXP pisSEXP, SEXP classesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type zs(zsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Ps(PsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pis(pisSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type classes(classesSEXP);
    rcpp_result_gen = Rcpp::wrap(update_z_single(zs, A, Ps, pis, classes));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mlsbm_logf_z_given_A", (DL_FUNC) &_mlsbm_logf_z_given_A, 4},
    {"_mlsbm_logf_z_given_A_single", (DL_FUNC) &_mlsbm_logf_z_given_A_single, 4},
    {"_mlsbm_rdirichlet_cpp", (DL_FUNC) &_mlsbm_rdirichlet_cpp, 2},
    {"_mlsbm_sample_SBM_fast", (DL_FUNC) &_mlsbm_sample_SBM_fast, 2},
    {"_mlsbm_update_P", (DL_FUNC) &_mlsbm_update_P, 5},
    {"_mlsbm_update_P_single", (DL_FUNC) &_mlsbm_update_P_single, 5},
    {"_mlsbm_update_counts", (DL_FUNC) &_mlsbm_update_counts, 2},
    {"_mlsbm_update_z", (DL_FUNC) &_mlsbm_update_z, 5},
    {"_mlsbm_update_z_single", (DL_FUNC) &_mlsbm_update_z_single, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_mlsbm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}