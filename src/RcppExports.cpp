// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// timesTwo
int timesTwo(int x);
RcppExport SEXP _NeverConsumersR_timesTwo(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(timesTwo(x));
    return rcpp_result_gen;
END_RCPP
}
// mod
/*  * Extend division reminder to vectors  *  * @param   a       Dividend   * @param   n       Divisor  */ double mod(double a, int n);
RcppExport SEXP _NeverConsumersR_mod(SEXP aSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(mod(a, n));
    return rcpp_result_gen;
END_RCPP
}
// gen_truncated_normals_never_c
arma::mat gen_truncated_normals_never_c(const arma::mat& trunc_value, const arma::mat& startxi, double numgen);
RcppExport SEXP _NeverConsumersR_gen_truncated_normals_never_c(SEXP trunc_valueSEXP, SEXP startxiSEXP, SEXP numgenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type trunc_value(trunc_valueSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type startxi(startxiSEXP);
    Rcpp::traits::input_parameter< double >::type numgen(numgenSEXP);
    rcpp_result_gen = Rcpp::wrap(gen_truncated_normals_never_c(trunc_value, startxi, numgen));
    return rcpp_result_gen;
END_RCPP
}
// update_Ni_with_covariates_c
Rcpp::List update_Ni_with_covariates_c(const arma::cube& Xtildei, const arma::mat& beta, const arma::mat& Utildei, const arma::vec& alpha, const arma::mat& GGalpha, const double n, const double mmi, const arma::vec& didconsume);
RcppExport SEXP _NeverConsumersR_update_Ni_with_covariates_c(SEXP XtildeiSEXP, SEXP betaSEXP, SEXP UtildeiSEXP, SEXP alphaSEXP, SEXP GGalphaSEXP, SEXP nSEXP, SEXP mmiSEXP, SEXP didconsumeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type Xtildei(XtildeiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Utildei(UtildeiSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type GGalpha(GGalphaSEXP);
    Rcpp::traits::input_parameter< const double >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double >::type mmi(mmiSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type didconsume(didconsumeSEXP);
    rcpp_result_gen = Rcpp::wrap(update_Ni_with_covariates_c(Xtildei, beta, Utildei, alpha, GGalpha, n, mmi, didconsume));
    return rcpp_result_gen;
END_RCPP
}
// gen_Wtildei_1foodplusenergy_never_c
arma::cube gen_Wtildei_1foodplusenergy_never_c(const arma::cube& Wtildei, const arma::mat& beta, const arma::cube& Xtildei, const arma::mat& Utildei, const double n, const arma::mat& iSigmae, const arma::mat& Wistar, const double mmi, const double numgen);
RcppExport SEXP _NeverConsumersR_gen_Wtildei_1foodplusenergy_never_c(SEXP WtildeiSEXP, SEXP betaSEXP, SEXP XtildeiSEXP, SEXP UtildeiSEXP, SEXP nSEXP, SEXP iSigmaeSEXP, SEXP WistarSEXP, SEXP mmiSEXP, SEXP numgenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type Wtildei(WtildeiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Xtildei(XtildeiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Utildei(UtildeiSEXP);
    Rcpp::traits::input_parameter< const double >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type iSigmae(iSigmaeSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Wistar(WistarSEXP);
    Rcpp::traits::input_parameter< const double >::type mmi(mmiSEXP);
    Rcpp::traits::input_parameter< const double >::type numgen(numgenSEXP);
    rcpp_result_gen = Rcpp::wrap(gen_Wtildei_1foodplusenergy_never_c(Wtildei, beta, Xtildei, Utildei, n, iSigmae, Wistar, mmi, numgen));
    return rcpp_result_gen;
END_RCPP
}
// formGofSigmae_never_c
double formGofSigmae_never_c(const double r, const double theta, const double s22, const double s33, const arma::cube& qq, const double mmi);
RcppExport SEXP _NeverConsumersR_formGofSigmae_never_c(SEXP rSEXP, SEXP thetaSEXP, SEXP s22SEXP, SEXP s33SEXP, SEXP qqSEXP, SEXP mmiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type s22(s22SEXP);
    Rcpp::traits::input_parameter< const double >::type s33(s33SEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type qq(qqSEXP);
    Rcpp::traits::input_parameter< const double >::type mmi(mmiSEXP);
    rcpp_result_gen = Rcpp::wrap(formGofSigmae_never_c(r, theta, s22, s33, qq, mmi));
    return rcpp_result_gen;
END_RCPP
}
// updated_parameter_r_never_c
double updated_parameter_r_never_c(const double r, const double theta, const double s22, const double s33, const arma::cube& qq, const double mmi, double n);
RcppExport SEXP _NeverConsumersR_updated_parameter_r_never_c(SEXP rSEXP, SEXP thetaSEXP, SEXP s22SEXP, SEXP s33SEXP, SEXP qqSEXP, SEXP mmiSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type s22(s22SEXP);
    Rcpp::traits::input_parameter< const double >::type s33(s33SEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type qq(qqSEXP);
    Rcpp::traits::input_parameter< const double >::type mmi(mmiSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(updated_parameter_r_never_c(r, theta, s22, s33, qq, mmi, n));
    return rcpp_result_gen;
END_RCPP
}
// updated_parameter_s22_never_c
double updated_parameter_s22_never_c(const double r, const double theta, const double s22, const double s33, const arma::cube& qq, const double mmi, double n);
RcppExport SEXP _NeverConsumersR_updated_parameter_s22_never_c(SEXP rSEXP, SEXP thetaSEXP, SEXP s22SEXP, SEXP s33SEXP, SEXP qqSEXP, SEXP mmiSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type s22(s22SEXP);
    Rcpp::traits::input_parameter< const double >::type s33(s33SEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type qq(qqSEXP);
    Rcpp::traits::input_parameter< const double >::type mmi(mmiSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(updated_parameter_s22_never_c(r, theta, s22, s33, qq, mmi, n));
    return rcpp_result_gen;
END_RCPP
}
// updated_parameter_s33_never_c
double updated_parameter_s33_never_c(const double r, const double theta, const double s22, const double s33, const arma::cube& qq, const double mmi, double n);
RcppExport SEXP _NeverConsumersR_updated_parameter_s33_never_c(SEXP rSEXP, SEXP thetaSEXP, SEXP s22SEXP, SEXP s33SEXP, SEXP qqSEXP, SEXP mmiSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type s22(s22SEXP);
    Rcpp::traits::input_parameter< const double >::type s33(s33SEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type qq(qqSEXP);
    Rcpp::traits::input_parameter< const double >::type mmi(mmiSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(updated_parameter_s33_never_c(r, theta, s22, s33, qq, mmi, n));
    return rcpp_result_gen;
END_RCPP
}
// updated_parameter_theta_never_c
double updated_parameter_theta_never_c(const double r, const double theta, const double s22, const double s33, const arma::cube& qq, const double mmi);
RcppExport SEXP _NeverConsumersR_updated_parameter_theta_never_c(SEXP rSEXP, SEXP thetaSEXP, SEXP s22SEXP, SEXP s33SEXP, SEXP qqSEXP, SEXP mmiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type s22(s22SEXP);
    Rcpp::traits::input_parameter< const double >::type s33(s33SEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type qq(qqSEXP);
    Rcpp::traits::input_parameter< const double >::type mmi(mmiSEXP);
    rcpp_result_gen = Rcpp::wrap(updated_parameter_theta_never_c(r, theta, s22, s33, qq, mmi));
    return rcpp_result_gen;
END_RCPP
}
// update_iSigmau_c
Rcpp::List update_iSigmau_c(const arma::mat& Sigmau, const double prior_Sigmau_doff, const arma::mat& prior_Sigmau_mean, const arma::mat& Utildei, const double n, double jjMCMC);
RcppExport SEXP _NeverConsumersR_update_iSigmau_c(SEXP SigmauSEXP, SEXP prior_Sigmau_doffSEXP, SEXP prior_Sigmau_meanSEXP, SEXP UtildeiSEXP, SEXP nSEXP, SEXP jjMCMCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Sigmau(SigmauSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_Sigmau_doff(prior_Sigmau_doffSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type prior_Sigmau_mean(prior_Sigmau_meanSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Utildei(UtildeiSEXP);
    Rcpp::traits::input_parameter< const double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type jjMCMC(jjMCMCSEXP);
    rcpp_result_gen = Rcpp::wrap(update_iSigmau_c(Sigmau, prior_Sigmau_doff, prior_Sigmau_mean, Utildei, n, jjMCMC));
    return rcpp_result_gen;
END_RCPP
}
// update_Utildei_c
arma::mat update_Utildei_c(arma::mat& Utildei, const arma::mat& beta, const arma::cube& Wtildei, const arma::mat& iSigmae, const arma::vec& Ni, const arma::vec& isnever, const arma::vec& didconsume, const arma::cube& Xtildei, const double mmi, const arma::mat& iSigmau, double n);
RcppExport SEXP _NeverConsumersR_update_Utildei_c(SEXP UtildeiSEXP, SEXP betaSEXP, SEXP WtildeiSEXP, SEXP iSigmaeSEXP, SEXP NiSEXP, SEXP isneverSEXP, SEXP didconsumeSEXP, SEXP XtildeiSEXP, SEXP mmiSEXP, SEXP iSigmauSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type Utildei(UtildeiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Wtildei(WtildeiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type iSigmae(iSigmaeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Ni(NiSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type isnever(isneverSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type didconsume(didconsumeSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Xtildei(XtildeiSEXP);
    Rcpp::traits::input_parameter< const double >::type mmi(mmiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type iSigmau(iSigmauSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(update_Utildei_c(Utildei, beta, Wtildei, iSigmae, Ni, isnever, didconsume, Xtildei, mmi, iSigmau, n));
    return rcpp_result_gen;
END_RCPP
}
// update_beta1_with_prior_mean_random_walk_c
arma::colvec update_beta1_with_prior_mean_random_walk_c(const arma::cube& Xtildei, const double mmi, const arma::mat& prior_beta_mean, const arma::cube prior_beta_cov, const arma::mat& beta, const arma::cube& Wtildei, const arma::mat& Utildei, const arma::mat& iSigmae, const arma::vec& isnever, const double update_beta1_var_ind);
RcppExport SEXP _NeverConsumersR_update_beta1_with_prior_mean_random_walk_c(SEXP XtildeiSEXP, SEXP mmiSEXP, SEXP prior_beta_meanSEXP, SEXP prior_beta_covSEXP, SEXP betaSEXP, SEXP WtildeiSEXP, SEXP UtildeiSEXP, SEXP iSigmaeSEXP, SEXP isneverSEXP, SEXP update_beta1_var_indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type Xtildei(XtildeiSEXP);
    Rcpp::traits::input_parameter< const double >::type mmi(mmiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type prior_beta_mean(prior_beta_meanSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type prior_beta_cov(prior_beta_covSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Wtildei(WtildeiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Utildei(UtildeiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type iSigmae(iSigmaeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type isnever(isneverSEXP);
    Rcpp::traits::input_parameter< const double >::type update_beta1_var_ind(update_beta1_var_indSEXP);
    rcpp_result_gen = Rcpp::wrap(update_beta1_with_prior_mean_random_walk_c(Xtildei, mmi, prior_beta_mean, prior_beta_cov, beta, Wtildei, Utildei, iSigmae, isnever, update_beta1_var_ind));
    return rcpp_result_gen;
END_RCPP
}
// update_beta1_with_prior_mean_c
arma::colvec update_beta1_with_prior_mean_c(const arma::cube& Xtildei, const double mmi, const arma::mat& prior_beta_mean, const arma::cube prior_beta_cov, const arma::mat& beta, const arma::cube& Wtildei, const arma::mat& Utildei, const arma::mat& iSigmae, const arma::colvec& isnever, const double update_beta1_var_ind);
RcppExport SEXP _NeverConsumersR_update_beta1_with_prior_mean_c(SEXP XtildeiSEXP, SEXP mmiSEXP, SEXP prior_beta_meanSEXP, SEXP prior_beta_covSEXP, SEXP betaSEXP, SEXP WtildeiSEXP, SEXP UtildeiSEXP, SEXP iSigmaeSEXP, SEXP isneverSEXP, SEXP update_beta1_var_indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type Xtildei(XtildeiSEXP);
    Rcpp::traits::input_parameter< const double >::type mmi(mmiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type prior_beta_mean(prior_beta_meanSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type prior_beta_cov(prior_beta_covSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Wtildei(WtildeiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Utildei(UtildeiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type iSigmae(iSigmaeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type isnever(isneverSEXP);
    Rcpp::traits::input_parameter< const double >::type update_beta1_var_ind(update_beta1_var_indSEXP);
    rcpp_result_gen = Rcpp::wrap(update_beta1_with_prior_mean_c(Xtildei, mmi, prior_beta_mean, prior_beta_cov, beta, Wtildei, Utildei, iSigmae, isnever, update_beta1_var_ind));
    return rcpp_result_gen;
END_RCPP
}
// update_beta2_with_prior_mean_c
arma::colvec update_beta2_with_prior_mean_c(const arma::cube& Xtildei, const double mmi, const arma::mat& prior_beta_mean, const arma::cube prior_beta_cov, const arma::mat& beta, const arma::cube& Wtildei, const arma::mat& Utildei, const arma::mat& iSigmae);
RcppExport SEXP _NeverConsumersR_update_beta2_with_prior_mean_c(SEXP XtildeiSEXP, SEXP mmiSEXP, SEXP prior_beta_meanSEXP, SEXP prior_beta_covSEXP, SEXP betaSEXP, SEXP WtildeiSEXP, SEXP UtildeiSEXP, SEXP iSigmaeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type Xtildei(XtildeiSEXP);
    Rcpp::traits::input_parameter< const double >::type mmi(mmiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type prior_beta_mean(prior_beta_meanSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type prior_beta_cov(prior_beta_covSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Wtildei(WtildeiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Utildei(UtildeiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type iSigmae(iSigmaeSEXP);
    rcpp_result_gen = Rcpp::wrap(update_beta2_with_prior_mean_c(Xtildei, mmi, prior_beta_mean, prior_beta_cov, beta, Wtildei, Utildei, iSigmae));
    return rcpp_result_gen;
END_RCPP
}
// update_beta3_with_prior_mean_c
arma::colvec update_beta3_with_prior_mean_c(const arma::cube& Xtildei, const double mmi, const arma::mat& prior_beta_mean, const arma::cube prior_beta_cov, const arma::mat& beta, const arma::cube& Wtildei, const arma::mat& Utildei, const arma::mat& iSigmae);
RcppExport SEXP _NeverConsumersR_update_beta3_with_prior_mean_c(SEXP XtildeiSEXP, SEXP mmiSEXP, SEXP prior_beta_meanSEXP, SEXP prior_beta_covSEXP, SEXP betaSEXP, SEXP WtildeiSEXP, SEXP UtildeiSEXP, SEXP iSigmaeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type Xtildei(XtildeiSEXP);
    Rcpp::traits::input_parameter< const double >::type mmi(mmiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type prior_beta_mean(prior_beta_meanSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type prior_beta_cov(prior_beta_covSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Wtildei(WtildeiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Utildei(UtildeiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type iSigmae(iSigmaeSEXP);
    rcpp_result_gen = Rcpp::wrap(update_beta3_with_prior_mean_c(Xtildei, mmi, prior_beta_mean, prior_beta_cov, beta, Wtildei, Utildei, iSigmae));
    return rcpp_result_gen;
END_RCPP
}
// ginverse_c
arma::mat ginverse_c(const arma::mat& z, const double lambda);
RcppExport SEXP _NeverConsumersR_ginverse_c(SEXP zSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(ginverse_c(z, lambda));
    return rcpp_result_gen;
END_RCPP
}
// backtransform_c
arma::mat backtransform_c(const double lambda, arma::mat& Xtildei, arma::vec& beta, double sigmae, double mumu, double sigsig, arma::vec& Utildei, double& n);
RcppExport SEXP _NeverConsumersR_backtransform_c(SEXP lambdaSEXP, SEXP XtildeiSEXP, SEXP betaSEXP, SEXP sigmaeSEXP, SEXP mumuSEXP, SEXP sigsigSEXP, SEXP UtildeiSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Xtildei(XtildeiSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type sigmae(sigmaeSEXP);
    Rcpp::traits::input_parameter< double >::type mumu(mumuSEXP);
    Rcpp::traits::input_parameter< double >::type sigsig(sigsigSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Utildei(UtildeiSEXP);
    Rcpp::traits::input_parameter< double& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(backtransform_c(lambda, Xtildei, beta, sigmae, mumu, sigsig, Utildei, n));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _NeverConsumersR_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _NeverConsumersR_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _NeverConsumersR_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _NeverConsumersR_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NeverConsumersR_timesTwo", (DL_FUNC) &_NeverConsumersR_timesTwo, 1},
    {"_NeverConsumersR_mod", (DL_FUNC) &_NeverConsumersR_mod, 2},
    {"_NeverConsumersR_gen_truncated_normals_never_c", (DL_FUNC) &_NeverConsumersR_gen_truncated_normals_never_c, 3},
    {"_NeverConsumersR_update_Ni_with_covariates_c", (DL_FUNC) &_NeverConsumersR_update_Ni_with_covariates_c, 8},
    {"_NeverConsumersR_gen_Wtildei_1foodplusenergy_never_c", (DL_FUNC) &_NeverConsumersR_gen_Wtildei_1foodplusenergy_never_c, 9},
    {"_NeverConsumersR_formGofSigmae_never_c", (DL_FUNC) &_NeverConsumersR_formGofSigmae_never_c, 6},
    {"_NeverConsumersR_updated_parameter_r_never_c", (DL_FUNC) &_NeverConsumersR_updated_parameter_r_never_c, 7},
    {"_NeverConsumersR_updated_parameter_s22_never_c", (DL_FUNC) &_NeverConsumersR_updated_parameter_s22_never_c, 7},
    {"_NeverConsumersR_updated_parameter_s33_never_c", (DL_FUNC) &_NeverConsumersR_updated_parameter_s33_never_c, 7},
    {"_NeverConsumersR_updated_parameter_theta_never_c", (DL_FUNC) &_NeverConsumersR_updated_parameter_theta_never_c, 6},
    {"_NeverConsumersR_update_iSigmau_c", (DL_FUNC) &_NeverConsumersR_update_iSigmau_c, 6},
    {"_NeverConsumersR_update_Utildei_c", (DL_FUNC) &_NeverConsumersR_update_Utildei_c, 11},
    {"_NeverConsumersR_update_beta1_with_prior_mean_random_walk_c", (DL_FUNC) &_NeverConsumersR_update_beta1_with_prior_mean_random_walk_c, 10},
    {"_NeverConsumersR_update_beta1_with_prior_mean_c", (DL_FUNC) &_NeverConsumersR_update_beta1_with_prior_mean_c, 10},
    {"_NeverConsumersR_update_beta2_with_prior_mean_c", (DL_FUNC) &_NeverConsumersR_update_beta2_with_prior_mean_c, 8},
    {"_NeverConsumersR_update_beta3_with_prior_mean_c", (DL_FUNC) &_NeverConsumersR_update_beta3_with_prior_mean_c, 8},
    {"_NeverConsumersR_ginverse_c", (DL_FUNC) &_NeverConsumersR_ginverse_c, 2},
    {"_NeverConsumersR_backtransform_c", (DL_FUNC) &_NeverConsumersR_backtransform_c, 8},
    {"_NeverConsumersR_rcpparma_hello_world", (DL_FUNC) &_NeverConsumersR_rcpparma_hello_world, 0},
    {"_NeverConsumersR_rcpparma_outerproduct", (DL_FUNC) &_NeverConsumersR_rcpparma_outerproduct, 1},
    {"_NeverConsumersR_rcpparma_innerproduct", (DL_FUNC) &_NeverConsumersR_rcpparma_innerproduct, 1},
    {"_NeverConsumersR_rcpparma_bothproducts", (DL_FUNC) &_NeverConsumersR_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_NeverConsumersR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
