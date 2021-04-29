// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include<cmath>
using namespace Rcpp;

// [[Rcpp::export]]
int timesTwo(int x)
{
  return x*2;
}

/***R
timesTwo(42)
*/

// [[Rcpp::export]]
/*
 * Extend division reminder to vectors
 *
 * @param   a       Dividend 
 * @param   n       Divisor
 */
double mod(double a, int n)
{
  return a - std::floor(a/n)*n;
}   


// [[Rcpp::export]]
arma::mat gen_truncated_normals_never_c(const arma::mat& trunc_value,
                                        const arma::mat& startxi,double numgen){
  
  // calling exp()
  //Function f("exp");   
  
  // arma::mat a(2,2);
  //   a.zeros();
  
  double  eps = 2.2204*std::exp(-16);
  int  n = trunc_value.n_rows;
  arma::mat alpha(n,1), aa(n,1), thesign(n,1);
  arma::mat genww = arma::zeros(n,1);
  for(int i = 0; i < n; i++) {
    alpha(i,0)  = (trunc_value(i,0) + std::sqrt(4 + std::pow(trunc_value(i,0), 2))) / 2;
    }
  //arma::mat alpha  = (trunc_value + std::sqrt(4 + std::pow(trunc_value, 2))) / 2;
  
  for(int i = 0; i< n; i++){
    if(trunc_value(i,0) >= 0.0)
      thesign(i,0) = 1.0;
    else
      thesign(i,0) = 0.0;
    }
  
  // Implementing arma::mat genww  = trunc_value % aa;
  
  for(int i = 0; i< n; i++){
    if(trunc_value(i,0) > 0.0)
      genww(i,0) = trunc_value(i,0);
  }
  
  arma::mat temp2  = arma::randn(n,1);
  arma::mat onemat = arma::ones(n,1);
  arma::mat mmmm(n,1), kkkk(n,1), hhhh(n,1);
  
  for(int jj = 1; jj<= numgen; jj++){
    arma::mat alpha_inv = arma::zeros(n,1);
    for (int i=0;i<n;i++){
      alpha_inv(i,0) = 1/alpha(i,0);
      }
    
    arma::mat unif = arma::randu(n,1);
    
    arma::mat term1 = alpha_inv % arma::log(unif);
   
    arma::mat xicand = trunc_value - term1; //( (1 / alpha) % log(arma::randu(n,1)));
    arma::mat gg = arma::exp(-0.5 * arma::pow( (xicand - alpha), 2));
    arma::mat ss = arma::randu(n,1);
    
    for(int i=0; i<n; i++){
      if( ss(i,0) < gg(i,0))
        mmmm(i,0) = 1.0;
      else
        mmmm(i,0) = 0.0;
      }

    
    arma::mat temp1  = (xicand % mmmm) + (genww % (onemat - mmmm));
    arma::mat ssss   = arma::randn(n,1);
    
    for(int i=0; i<n; i++){
      if( ssss(i,0) < trunc_value(i,0)){
        kkkk(i,0) = 1.0;
      }else{
        kkkk(i,0) = 0.0;
      }}
    temp2  = (temp2 % kkkk) +       (ssss % (onemat - kkkk));
    genww  = (temp2 % (onemat - thesign) + (temp1 % thesign ));
  }
 
  for(int i=0; i<n; i++){
    if( genww(i,0) > trunc_value(i,0)){
      hhhh(i,0) = 1.0;
    }else{
      hhhh(i,0) = 0.0;
    }
  }
  arma::mat eps_mat = arma::zeros(n,1);
  for(int i=0; i<n; i++){
    eps_mat(i,0) = eps;
    }
  
  genww  = (genww % hhhh) + ((trunc_value + eps_mat) % (onemat - hhhh));
  return(genww);
//return(arma::ones(10));
  
}


// [[Rcpp::export]]
Rcpp::List update_Ni_with_covariates_c(const arma::cube& Xtildei,const arma::mat& beta,
                                       const arma::mat& Utildei,const arma::vec& alpha,
                                       const arma::mat& GGalpha,const double n,
                                       const double mmi,const arma::vec& didconsume){

  arma::mat  tt = (Xtildei.slice(0) * beta.col(0)) + Utildei.col(0);
  //std::cout << "Xtildei.slice(1)" << Xtildei.slice(1) <<std::endl;
  //std::cout << "beta.col(1)" << beta.col(1) <<std::endl;
  //std::cout << "Xtildei.slice(1)" << Xtildei.slice(1) <<std::endl;
  //std::cout << "tt" << tt <<std::endl;
  arma::mat aq1 = arma::normcdf(tt,0,1);
  arma::mat aq2 = arma::normcdf(GGalpha*alpha,0,1);
  //std::cout << "aq1" << aq1 <<std::endl;
  //std::cout << "aq2" << aq2 <<std::endl;
  arma::mat denom = arma::zeros(aq2.n_rows,1);
  arma::mat midterm(n,1);
  for(int i=0;i<aq2.n_rows;i++){
    denom(i,0) = 1/(1-aq2(i,0));
    midterm(i,0) = std::pow((1 - aq1(i,0)), mmi);
    }
  arma::mat cc1 = aq2 % midterm % denom;
  //std::cout << "CC1" << cc1 <<std::endl;
  arma::mat onemat = arma::ones(n,1);
  arma::mat ppi(cc1.n_rows, 1) ;
  for(int i = 0; i<cc1.n_rows; i++){
    ppi(i,0) =  1 / (1 + cc1(i, 0));
    }
  arma::mat Ni = arma::zeros(n,1);
  arma::mat kk = GGalpha * alpha;
  //std::cout << "kk" << kk <<std::endl;
  arma::mat genww1  = gen_truncated_normals_never_c(-kk,-kk % arma::ones(n,1),50.0);
  arma::mat genww2  = gen_truncated_normals_never_c(kk,-kk % arma::ones(n,1),50.0);

  arma::mat uu = arma::randu(n,1);
  arma::mat rri(n,1);
  for(int i=0; i< n; i++){
    if(uu(i,0) < ppi(i,0)){
      rri(i,0) = 1.0;
    }else{
      rri(i,0) = 0.0;
    }
  }


  Ni = GGalpha*alpha + (didconsume % genww1) +
    ((onemat-didconsume) % (((onemat-rri) % genww1) - (rri % genww2)));

  return Rcpp::List::create( Rcpp::Named("Ni") = Ni, Rcpp::Named("ppi") = ppi);
}

// [[Rcpp::export]]
arma::cube gen_Wtildei_1foodplusenergy_never_c(const arma::cube& Wtildei,const arma::mat& beta,const arma::cube& Xtildei,
                                               const arma::mat& Utildei,const double n,const arma::mat& iSigmae,const arma::mat& Wistar,const double mmi,const double numgen){

  arma::cube Wtildeinew  = Wtildei;
  arma::cube Wtildeio    = Wtildei;
  double C2;
  arma::vec C1(beta.n_rows), qq(beta.n_rows);
  int varnum      = 0;
  for(int kk = 0; kk< mmi; kk++){
    C2      = 1 / iSigmae(varnum,varnum);
    C1      = iSigmae(varnum,varnum) * ((Xtildei.slice(varnum) * beta.col(varnum)) + Utildei.col(varnum));
    // std::cout<< "y is"<< ((Xtildei.slice(varnum) * beta.col(varnum)) + Utildei.col(varnum)) <<std::endl;
    // std::cout<< "iSigmae(varnum,varnum)" << iSigmae(varnum,varnum) << std::endl;
    // std::cout << C1 <<std::endl;
    // std::cout << (C2*C1) << std::endl;
    for(int jj = 0; jj< Xtildei.n_slices; jj++){
      if(abs(jj - varnum) > 0) {
        arma::mat w = Wtildei.slice(kk);

        qq = (w.col(jj) - (Xtildei.slice(jj) * beta.col(jj)) - Utildei.col(jj));
        // std::cout << "C1 initial" << C1 <<std::endl;
        // std::cout << "qq" << qq <<std::endl;
        C1 = C1 - (iSigmae(varnum,jj) * qq);
        // std::cout << "iSigmae(varnum,jj)" << iSigmae(varnum,jj) <<std::endl;
        // std::cout << "C1 new" << C1 <<std::endl;
      }

    }
    arma::vec mu      = C2 * C1;
    double sigma   = sqrt(C2);
    arma::vec startxi = mu/sigma;
    // std::cout << "mu" << mu <<std::endl;
    // std::cout << "sigma" << sigma <<std::endl;
    // std::cout << "startxi" << startxi <<std::endl;
    arma::vec genww1  = gen_truncated_normals_never_c(-mu/sigma,-startxi,numgen);
    arma::vec genww2  = gen_truncated_normals_never_c(mu/sigma,-startxi,numgen);
    // std::cout << "genww1" << genww1 <<std::endl;
    // std::cout << "genww2" << genww2 <<std::endl;
    //arma::mat w = Wtildeinew.slice(kk);
    Wtildeinew(arma::span::all,arma::span(varnum),arma::span(kk)) = mu + (sigma * ((Wistar.col(kk) % genww1) - ((1 - Wistar.col(kk)) % genww2)));

  }


  varnum      = 1;
  for(int kk = 0; kk< mmi; kk++){
    C2      = 1 / iSigmae(varnum,varnum);
    C1      = iSigmae(varnum,varnum) * ((Xtildei.slice(varnum) * beta.col(varnum)) + Utildei.col(varnum));
    for(int jj = 0; jj< Xtildei.n_slices; jj++){
      if(abs(jj - varnum) > 0) {
        arma::mat w = Wtildei.slice(kk);

        qq = (w.col(jj) - (Xtildei.slice(jj) * beta.col(jj)) - Utildei.col(jj));
        // std::cout << "C1 initial" << C1 <<std::endl;
        // std::cout << "qq" << qq <<std::endl;
        C1 = C1 - (iSigmae(varnum,jj) * qq);
        // std::cout << "iSigmae(varnum,jj)" << iSigmae(varnum,jj) <<std::endl;
        // std::cout << "C1 new" << C1 <<std::endl;
      }

    }
    arma::vec mu      = C2 * C1;
    double sigma   = sqrt(C2);
    // std::cout << "mu" << mu <<std::endl;
    // std::cout << "sigma" << sigma <<std::endl;
    arma::mat check = arma::randn(n,1);
    // std::cout << "X" << check <<std::endl;
    // std::cout << "mu + sigma*X" << mu + (sigma * check) <<std::endl;
    //arma::mat w = Wtildeinew.slice(kk);
    Wtildeinew(arma::span::all,arma::span(varnum),arma::span(kk)) = mu + (sigma * arma::randn(n,1));
  }

  arma::mat w = Wtildei.col(varnum);
  arma::mat wn = Wtildeinew.col(varnum);
  
  Wtildeinew.col(varnum) = (w % Wistar) + (wn % (1-Wistar));

  return(Wtildeinew) ;
}

// [[Rcpp::export]]
double formGofSigmae_never_c(const double r, const double theta, const double s22, const double s33, const arma::cube& qq,const double mmi){

  // Compute the loglikelihood contribution for Sigmae other than its
  // determinant
  arma::mat R = {{1, 0, r*cos(theta)},
  {0,  1, r*sin(theta)},
  {r*cos(theta), r*sin(theta), 1}};
  arma::vec v = {1, sqrt(s22), sqrt(s33)};
  arma::mat A = arma::diagmat(v);

  arma::mat Sigmae  = A * R * A;
  //std::cout << "C function Sigmae" << Sigmae <<std::endl;
  arma::mat iSigmae = inv(Sigmae);

  arma::mat tempMat(qq.n_cols,qq.n_cols);
  tempMat.zeros();

  for (int i=0; i< mmi; i++){
    tempMat = tempMat + arma::trans(qq.slice(i)) * qq.slice(i);
  }
  tempMat = iSigmae % tempMat;
  double  GofSigmae = -0.5 * (arma::accu(tempMat));

  return(GofSigmae);
}

// [[Rcpp::export]]
double updated_parameter_r_never_c(const double r, const double theta, const double s22, const double s33, const arma::cube& qq,const double mmi, double n){

  // Do the Metropolis Step for r.
  // r = current value
  // theta = current value
  // s22 = current value
  // s33 = current value
  double rcurr, rcand, ss, ss1, rnew, a, b, c, d, i;

  rcurr = r;
  arma::vec  rpossible = arma::linspace(-0.99, 0.99, 41);
  double  spacing   = rpossible(2) - rpossible(1);
  // if(abs(round(r,2)) <= 0.99 ){
  //   rcurr = round(r,2);
  // }else if(round(r,2) > 0.99){
  //   rcurr = 0.99;
  // }else{
  //   rcurr = -0.99;
  // }

  ss    = arma::randu<double>();
  if(ss <= 0.33){
    a = 1;
  }else{
    a = 0;
  }
  if(ss > 0.33 && ss <= 0.66){
    b = 1;
  }
  else{
    b = 0;
  }
  if(ss > 0.66){
    c = 1;
  }
  else{
    c = 0;
  }

  if (rcurr == -0.99){
    i = 1.0;
    rcand = (rcurr * a) + ((rcurr + spacing) *b) + ((rcurr + (2 * spacing)) *c);
   //std::cout << "this is ss:" << ss << " when rcurr = "<< rcurr << " and rcand = " << rcand << " and spacing = "<< spacing<< " and i = "<< i <<std::endl;

  }
  if (rcurr == 0.99){
   i = 2.0;
    rcand = (rcurr * a) + ((rcurr - spacing)*b) + ((rcurr - (2 * spacing)) * c);
    //std::cout << "this is ss:" << ss << " when rcurr = "<< rcurr << " and rcand = " << rcand << " and spacing = "<< spacing<< " and i = "<< i <<std::endl;
    
  }
  if (rcurr > -0.99 && rcurr < 0.99){
    i = 3.0;
    rcand = (rcurr * a) + ((rcurr + spacing)*b) + ((rcurr - spacing) * c);
    //std::cout << "this is ss:" << ss << " when rcurr = "<< rcurr << " and rcand = " << rcand << " and spacing = "<< spacing<< " and i = "<< i <<std::endl;
    
  }

  double GofSigmaecurr = formGofSigmae_never_c(rcurr,theta,s22,s33,qq,mmi);
  double GofSigmaecand = formGofSigmae_never_c(rcand,theta,s22,s33,qq,mmi);

  double gg            = GofSigmaecand - GofSigmaecurr;
  gg            = gg - ((mmi*n/2) * log(1 - pow(rcand, 2))) + ((mmi*n/2)
                                                                 * log(1 - pow(rcurr, 2)));
  gg            = std::min(1.0, exp(gg));
  ss1            = arma::randu<double>();

  if(ss1 < gg){
    d = 1;
  }else{
    d = 0;
    }
  rnew          = (rcand * d) + ((rcurr * (1 - d)));
  //std::cout << "this is ss1 = " << ss1 << " when gg = "<< gg << " then rnew = " << rnew <<std::endl;
  
  return(rnew);
}

// [[Rcpp::export]]
double updated_parameter_s22_never_c(const double r, const double theta, const double s22, const double s33, const arma::cube& qq,const double mmi, double n){

  // # Do the Metropolis Step for s22.
  // # r = current value
  // # theta = current value
  // # s22 = current value
  // # s33 = current value
  double a, d, ss, s22new;
  double s22curr       = s22;
  double  s22cand       = s22 + (0.4 * (arma::randu<double>() - 0.5));
  
  //std::cout<<"S22curr = "<< s22curr << " and s22cand = " <<s22cand<<std::endl;
  double GofSigmaecurr = formGofSigmae_never_c(r,theta,s22curr,s33,qq,mmi);
  double GofSigmaecand = formGofSigmae_never_c(r,theta,s22cand,s33,qq,mmi);

  double  gg            = GofSigmaecand - GofSigmaecurr;
  gg            = gg - ((mmi*n/2) * log(s22cand)) + ((mmi*n/2) * log(s22curr));
  //std::cout<<"gg initial = "<< gg <<std::endl;
  if(s22cand>=0 && s22cand<=3){
    a = 1;
  }else{
    a = 0;
  }
  gg            = std::min(1.0, exp(gg)*a);
  ss            = arma::randu<double>();
  // std::cout<<"gg final = "<< gg << " and ss = " <<ss <<std::endl;
  if(ss < gg){
    d = 1;
  }else{
    d = 0;
    }
  s22new       = (s22cand * d) +  (s22curr * (1 - d));
  // std::cout << "s22new = "<< s22new <<std::endl;

  return(s22new);

}

// [[Rcpp::export]]
double updated_parameter_s33_never_c(const double r, const double theta, const double s22, const double s33, const arma::cube& qq,const double mmi, double n){

  // # Do the Metropolis Step for s33.
  // # r = current value
  // # theta = current value
  // # s22 = current value
  // # s33 = current value

  double a, d, ss, s33new;
  double s33curr       = s33;
  double s33cand       = s33 + (0.4 * (arma::randu<double>() - 0.5));
  //std::cout<<"S33curr = "<< s33curr << " and s33cand = " <<s33cand<<std::endl;
  
  double GofSigmaecurr = formGofSigmae_never_c(r,theta,s22,s33curr,qq,mmi);
  double GofSigmaecand = formGofSigmae_never_c(r,theta,s22,s33cand,qq,mmi);

  double gg            = GofSigmaecand - GofSigmaecurr;
  gg            = gg - ((mmi*n/2) * log(s33cand)) + ((mmi*n/2) * log(s33curr));
  //std::cout<<"gg initial = "<< gg <<std::endl;
  
  if(s33cand>=0 && s33cand<=3){
    a = 1;
  }else{
    a = 0;
  }

  gg            = std::min(1.0, exp(gg)*a);

  ss            = arma::randu<double>();
  
  //std::cout<<"gg final = "<< gg << " and ss = " <<ss <<std::endl;
  if(ss < gg){
    d = 1;
  }else{
    d = 0;
    }
  s33new       = (s33cand * d) +  (s33curr * (1-d));
  //std::cout << "s33new = "<< s33new <<std::endl;
  return(abs(s33new));

}

// [[Rcpp::export]]
double updated_parameter_theta_never_c(const double r, const double theta, const double s22, const double s33, const arma::cube& qq,const double mmi){

  // # Do the Metropolis Step for theta.
  // # r = current value
  // # theta = current value
  // # s22 = current value
  // # s33 = current value
  double thetamin, thetamax,spacing, thetacurr, thetacand, a, b, c, d, ss;
  arma::vec thetapossible = 3.141593 * arma::linspace(-0.99, 0.99, 41);
  double theta_min = arma::min(thetapossible);
  double theta_max = arma::max(thetapossible);
  thetamin      = std::min(theta_min, theta);
  thetamax      = std::max(theta_max, theta);
  spacing       = thetapossible(2) - thetapossible(1);
  thetacurr     = theta;

  ss    = arma::randu<double>();
  if(ss <= 0.33){
    a = 1;
  }else{
    a = 0;
  }
  if(ss > 0.33 && ss <= 0.66){
    b = 1;
  }
  else{
    b = 0;
  }
  if(ss > 0.66){
    c = 1;
  }
  else{
    c = 0;
  }

  if(thetacurr <= thetamin){

    thetacand = (thetacurr * a) + ((thetacurr + spacing) *b )+
      ((thetacurr + (2 * spacing)) * c);
  }
  if(thetacurr >= thetamax){

    thetacand = (thetacurr * a) + ((thetacurr - spacing)*b)+
      ((thetacurr - (2 * spacing)) * c);
  }
  if (thetacurr  > thetamin){
    if (thetacurr < thetamax){

      thetacand = (thetacurr * a) + ((thetacurr + spacing) * b)
      + ((thetacurr - spacing) * c);
    }
  }
  double  GofSigmaecurr = formGofSigmae_never_c(r,thetacurr,s22,s33,qq,mmi);
  double  GofSigmaecand = formGofSigmae_never_c(r,thetacand,s22,s33,qq,mmi);

  double  gg            = std::min(1.0, exp(GofSigmaecand - GofSigmaecurr));
  ss            = arma::randu<double>();
  if(ss < gg){
    d = 1;
  }else{
    d = 0;
    }
  double thetanew      = (thetacand * d) + (thetacurr * (1-d));

  return(thetanew);
}

// [[Rcpp::export]]
Rcpp::List update_iSigmau_c(const arma::mat& Sigmau,const double prior_Sigmau_doff,
                            const arma::mat& prior_Sigmau_mean,const arma::mat& Utildei,
                            const double n, double jjMCMC){

  arma::mat aa  = ( (prior_Sigmau_doff - Sigmau.n_rows - 1) * prior_Sigmau_mean) +
    (Utildei.t() * Utildei);
  double      bb         = prior_Sigmau_doff + n;

  // std::cout<<aa<<std::endl;
  // std::cout<<bb<<std::endl;
  // 
  aa  = aa + (.00001 * (1 * n) * arma::eye(3, 3));
  bb  = bb + (.00001 * (1 * n));
  aa  = (aa + aa.t()) / (2);
  // std::cout<<aa<<std::endl;
  // std::cout<<bb<<std::endl;
  

  arma::vec eigval = arma::eig_sym(aa/bb);
  double thecount2 = 0.0;
  
  //std::cout<<aa<<std::endl;
  //std::cout<<eigval<<std::endl;
  
  arma::mat Sigmau_new(Sigmau.n_rows, Sigmau.n_cols), iSigmau_new(Sigmau.n_rows, Sigmau.n_cols);
  if( arma::min(eigval) <= 0.001){
    thecount2 = thecount2 + 1;
    if( thecount2 == 1){
      std::cout<<'Problem with Sigmau at step '<<jjMCMC<<std::endl;
    }
  }
  if(arma::min(eigval) > 0.001) {
    Sigmau_new     = arma::iwishrnd(aa,bb);
    iSigmau_new    = inv(Sigmau_new);

  }

  return Rcpp::List::create( Rcpp::Named("iSigmau_new") = iSigmau_new, Rcpp::Named("Sigmau_new") = Sigmau_new );
}

// [[Rcpp::export]]
arma::mat update_Utildei_c(arma::mat& Utildei,const arma::mat& beta,
                           const arma::cube& Wtildei,const arma::mat& iSigmae,
                           const arma::vec& Ni,const arma::vec& isnever,
                           const arma::vec& didconsume,const arma::cube& Xtildei,
                           const double mmi,const arma::mat& iSigmau,double n){

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //   % Update Utildei. This is done in two steps. In the first step, we generate
  //   % it assuming that everyone is a consumer. In the second step, those who
  //   % are never consumers, i.e., Ni < 0, have their values updated by a
  //   % Metropolis step.
  //   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //   %
  //   % INPUT
  //   % Utildei        = the latent person-specific effects
  //   % beta           = the current values of all the betas
  //   % Wtildei        = the latent and non-latent responses
  //   % iSigmae        = inverse of Sigmae
  //   % Ni             = the crucial random variable
  //   % isnever        = Indicator that Ni < 0 and the person is a never consumer
  //   % didconsume     = Indicator that the person consumed on one of the
  //   %                     recalls
  //   % Xtildei        = design matrix
  //   % mmi            = # of recalls (assumed same for all)
  //   % iSigmau        = inverse of Sigmaeu
  //
  //
  //   % Here is Step 1
  arma::mat  Utildei_Current = Utildei;

  arma::mat qq = (Xtildei.slice(0) * beta.col(0));
  for(int jj = 1; jj<=2; jj++){
    qq = arma::join_rows(qq, (Xtildei.slice(jj) * beta.col(jj)));
  }

  arma::cube y(qq.n_rows, qq.n_cols, mmi);
  for(int i=0; i< mmi; i++){
    y.slice(i) = qq;
  }
  arma::cube  ss = Wtildei - y;
  //std::cout << "ss" << ss << std::endl;
  
  arma::mat f = arma::sum(ss, 2);
  //std::cout << "f" << f << std::endl;
  
  arma::mat c1  = arma::trans(iSigmae * arma::trans(f));
  //std::cout << "c1" << c1 << std::endl;
  arma::mat c2  = inv(iSigmau + (mmi * iSigmae));
  // std::cout<<"iSigmau"<<iSigmau<<std::endl;
  // std::cout<<"iSigmae"<<iSigmae<<std::endl;
  // std::cout << "c2" << c2 << std::endl;
  // 
  // // Here is Step 2

  //c2 = (c2 + arma::trans(c2))/2;
  //std::cout<< "(c2+ c2')/2 = " << c2 << std::endl;
  arma::mat g = arma::sqrtmat_sympd(c2);
  arma::vec  one_vec = arma::ones(n);
  //std::cout << "g = " << g << std::endl;
  Utildei = (c2 * c1.t()).t() + (arma::randn(n,3) * g);
  arma::vec c2_cand = 1 / ( arma::pow(one_vec - arma::normcdf(qq.col(0) + Utildei.col(0)), mmi));
  arma::vec c2_curr = 1 / ( arma::pow(one_vec - arma::normcdf(qq.col(0) + Utildei_Current.col(0)), mmi));
  
  // std::cout << "1/c2_cand" << ( arma::pow(1 - arma::normcdf(qq.col(0) + Utildei.col(0)), mmi)) <<std::endl;
  // std::cout << "c2_cand" << c2_cand <<std::endl;
  
  arma::vec uu = arma::randu(n);
  int p = c2_cand.n_rows;
  arma::vec rri(p);
  for(int i=0; i< p; i++){
    if(uu[i] < (c2_cand[i] / c2_curr[i])){
      rri[i] = 1;
    }else{
      rri[i] = 0;
    }

  }

  arma::vec  one = arma::ones(p);
  double ddff       = Utildei.n_cols;
  arma::mat Utildei_MH = (Utildei % (rri * arma::ones(1,ddff)) ) +
                          (Utildei_Current % ((one - rri) * arma::ones(1,ddff)));

  arma::mat Utildei_new = (Utildei % ((one - isnever) * arma::ones(1,ddff))) + 
                          (Utildei_MH % (isnever * arma::ones(1,ddff)));

  return(Utildei_new);
}


// [[Rcpp::export]]
arma::colvec update_beta1_with_prior_mean_random_walk_c(const arma::cube& Xtildei, 
                                                        const double mmi,
                                                        const arma::mat& prior_beta_mean,
                                                        const arma::cube prior_beta_cov, 
                                                        const arma::mat& beta,
                                                        const arma::cube& Wtildei, 
                                                        const arma::mat& Utildei,
                                                        const arma::mat& iSigmae, 
                                                        const arma::vec& isnever, 
                                                        const double update_beta1_var_ind){

  arma::mat xx         = Xtildei.slice(0);
  arma::mat cc2        = inv(inv(prior_beta_cov.slice(0)) + (mmi * iSigmae(0,0) * (xx.t() * xx)));
  double mmbeta        = beta.n_rows;
  arma::mat cc1        = arma::zeros(mmbeta,1);
  cc1        = cc1 + (inv(prior_beta_cov.slice(0)) * prior_beta_mean.col(0));
  for( int jji=0; jji< mmi; jji ++){
    arma::mat w = Wtildei.slice(jji);
    cc1    = cc1 + (iSigmae(0,0)*(xx.t() * (w.col(0) - Utildei.col(0))));
    cc1    = cc1 + (iSigmae(0,1)*(xx.t() * (w.col(1) - (xx * beta.col(1)) - Utildei.col(1))));
    cc1    = cc1 + (iSigmae(0,2)*(xx.t() * (w.col(2) - (xx * beta.col(2)) - Utildei.col(2))));
  }

  arma::vec beta1_curr = beta.col(0);
  cc2 = (cc2 + cc2.t())/2;
  // std::cout<<"cc2 = "<< cc2 << std::endl;
  //std::cout<<"cc2/update_beta1_var_ind = " << cc2/update_beta1_var_ind<<std::endl;
  arma::mat t = arma::sqrtmat_sympd(cc2/update_beta1_var_ind);
  //std::cout<< "t = " << t << std::endl;
  arma::mat random = arma::randn(mmbeta,1);
  //std::cout<< "random = " << random << std::endl;
  //std::cout<< "t * random = " << t * random << std::endl;
  arma::vec beta1_cand = beta1_curr + (t * random);
  //std::cout<< "beta1_cand = " << beta1_cand << std::endl;
  //std::cout<< "beta1_curr = " << beta1_curr << std::endl;
  double n = xx.n_rows;
  arma::vec  v = arma::ones(n);
 
  double lc2_cand   = -mmi * sum((isnever % log((v - arma::normcdf((xx * beta1_cand) + Utildei.col(0))))));
  double lc2_curr   = -mmi * sum((isnever % log((v - arma::normcdf((xx * beta1_curr) + Utildei.col(0))))));
  // std::cout <<"cc1 size" << cc1.n_rows << "times" << cc1.n_cols << std::endl;
  // std::cout <<"beta1 size" << beta1_cand.n_rows << "times" << beta1_cand.n_cols << std::endl;
  arma::mat lc1_cand   = cc1.t() * beta1_cand - beta1_cand.t() * inv(cc2) * beta1_cand /2;
  arma::mat lc1_curr   = cc1.t() * beta1_curr - beta1_curr.t() * inv(cc2) * beta1_curr /2;
  //std::cout <<"Check 1" << std::endl;
  // std::cout << "lc1_cand" << lc1_cand << std::endl;
  // std::cout << "lc1_cand size = " << lc1_cand.n_rows<<"times"<<lc1_cand.n_cols << std::endl;
  arma::mat A(1,1);
  A.ones();
  arma::mat gghh       = arma::min(A ,exp( lc1_cand + lc2_cand - lc1_curr - lc2_curr));
  // std::cout << "lc1_cand = "<< lc1_cand << " lc2_cand = "<< lc2_cand << " lc1_curr = " << lc1_curr << " lc2_curr = "<<lc2_curr << std::endl;
  // std::cout << "lc1_cand + lc2_cand - lc1_curr - lc2_curr" << lc1_cand + lc2_cand - lc1_curr - lc2_curr << std::endl;
  // std::cout<< "exp(lc1_cand + lc2_cand - lc1_curr - lc2_curr) = " << exp(lc1_cand + lc2_cand - lc1_curr - lc2_curr)<<std::endl;
  // std::cout << "gghh" << gghh <<  std::endl;
  arma::mat ss(1,1);
  ss(0,0)            = arma::randu<double>();
  arma::vec rri(beta1_cand.n_elem);
  if(ss(0,0) < gghh(0,0)){
    rri = arma::ones(beta1_cand.n_elem);
  }else{
    rri = arma::zeros(beta1_cand.n_elem);}
  arma::vec  one = arma::ones(beta.n_rows);
  arma::vec beta1      = (beta1_cand % rri) + (beta1_curr % (one - rri));

  return(beta1);

}


// [[Rcpp::export]]
arma::colvec update_beta1_with_prior_mean_c(const arma::cube& Xtildei, const double mmi, const arma::mat& prior_beta_mean,
                                            const arma::cube prior_beta_cov, const arma::mat& beta, const arma::cube& Wtildei, const arma::mat& Utildei,
                                            const arma::mat& iSigmae, const arma::colvec& isnever, const double update_beta1_var_ind){

  // # This program updated beta1 using  Metropolis step. It also uses the fact
  // # that the design matrix for the food is the same for every repeat. This
  // # code would have to be changed is this were not true.
  // #
  // # Input
  // # Xtildei        = design matrix
  // # mmi            = # of recalls (assumed same for all)
  // # cov_prior_beta = the prior covariance matrix for beta1
  // # beta           = the current values of all the betas
  // # Wtildei        = the latent and non-latent responses
  // # Utildei        = the latent person-specific effects
  // # iSigmae        = the 2x2 inverse of Sigmae. iSigmae(1,1) = 1
  arma::mat xx         = Xtildei.slice(0);
  arma::mat cc2        = arma::inv(arma::inv(prior_beta_cov.slice(0)) +
    (mmi * iSigmae(0,0) * (xx.t()* xx)));
  double mmbeta     = beta.n_rows;
  arma::mat cc1        = arma::zeros(mmbeta,1);
  cc1        = cc1 + (arma::inv(prior_beta_cov.slice(0)) * prior_beta_mean.col(0));
  for (int jji = 0; jji< mmi; jji++){
    arma::mat w = Wtildei.slice(jji);
    cc1    = cc1 + (iSigmae(0,0)*(xx.t() * (w.col(0) - Utildei.col(0))));
    cc1    = cc1 + (iSigmae(0,1)*(xx.t() * (w.col(1) - (xx * beta.col(1)) - Utildei.col(1))));
    cc1    = cc1 + (iSigmae(0,2)*(xx.t() * (w.col(2) - (xx * beta.col(2)) - Utildei.col(2))));
  }
  //cc2 = (cc2 + cc2.t())/2;
  arma::mat t = arma::sqrtmat_sympd(cc2/update_beta1_var_ind);
  arma::vec beta1_cand = (cc2 * cc1) + (t* arma::randn(mmbeta,1));
  arma::vec  beta1_curr = beta.col(0);
  double n = xx.n_rows;
  arma::vec  v = arma::ones(n);
  double lc2_cand   = -mmi * sum((isnever %log((v - arma::normcdf((xx * beta1_cand) + Utildei.col(0))))));
  double lc2_curr   = -mmi * sum((isnever %log((v - arma::normcdf((xx* beta1_curr) + Utildei.col(0))))));
  arma::mat lc1_cand   = cc1.t() * beta1_cand - beta1_cand.t() * inv(cc2) * beta1_cand /2;
  arma::mat lc1_curr   = cc1.t() * beta1_curr - beta1_curr.t() * inv(cc2) * beta1_curr /2;

  arma::mat A(1,1);
  A.ones();

  arma::mat gghh       = arma::min(A ,exp( (1-update_beta1_var_ind) *lc1_cand + lc2_cand
                                             - (1-update_beta1_var_ind) * lc1_curr - lc2_curr));

  arma::mat ss(1,1);
  ss(0,0)            = arma::randu<double>();
  arma::vec rri(beta.n_rows);
  if(ss(0,0) < gghh(0,0)){
    rri = arma::ones(beta.n_rows);
  }else{
    rri = arma::zeros(beta.n_rows);}
  arma::vec  one = arma::ones(beta.n_rows);
  arma::colvec beta1      = (beta1_cand % rri) + (beta1_curr % (one - rri));

  return(beta1);

}

// [[Rcpp::export]]
arma::colvec update_beta2_with_prior_mean_c(const arma::cube& Xtildei, const double mmi, const arma::mat& prior_beta_mean,
                                            const arma::cube prior_beta_cov, const arma::mat& beta, const arma::cube& Wtildei, const arma::mat& Utildei,
                                            const arma::mat& iSigmae){

  // % This program updated beta2, and requires no Metropolis step. It also uses
  // % the fact that the design matrix for the food is the same for every
  // % repeat. This code would have to be changed is this were not true.
  // %
  // % There is another special feature about this problem, namely that since
  // % there is no energy, Sigmae is diagonal. This makes the step for beta2 a
  // % lot different.
  // %
  // % Input
  // % Xtildei        = design matrix
  // % mmi            = # of recalls (assumed same for all)
  // % cov_prior_beta = the prior covariance matrix for beta1
  // % beta           = the current values of all the betas
  // % Wtildei        = the latent and non-latent responses
  // % Utildei        = the latent person-specific effects
  // % iSigmae        = the 2x2 inverse of Sigmae. iSigmae(1,1) = 1

  arma::mat xx         = Xtildei.slice(0);
  arma::mat cc2        = inv(inv(prior_beta_cov.slice(1)) + (mmi * iSigmae(1,1) * (xx.t() * xx)));
  double mmbeta        = beta.n_rows;
  arma::mat cc1        = arma::zeros(mmbeta,1);
  cc1        = cc1 + (inv(prior_beta_cov.slice(1)) * prior_beta_mean.col(1));
  for( int jji=0; jji< mmi; jji ++){
    arma::mat w = Wtildei.slice(jji);
    cc1    = cc1 + (iSigmae(1,1)*(xx.t() * (w.col(1) - Utildei.col(1))));
    cc1    = cc1 + (iSigmae(0,1)*(xx.t() * (w.col(0) - (xx * beta.col(0)) - Utildei.col(0))));
    cc1    = cc1 + (iSigmae(1,2)*(xx.t() * (w.col(2) - (xx * beta.col(2)) - Utildei.col(2))));
  }
  //cc2 = (cc2 + cc2.t())/2;
  arma::mat t = arma::sqrtmat_sympd(cc2);
  arma::colvec beta2      = (cc2 * cc1) + (t * arma::randn(mmbeta,1));
  return(beta2);
}

// [[Rcpp::export]]
arma::colvec update_beta3_with_prior_mean_c(const arma::cube& Xtildei, const double mmi, const arma::mat& prior_beta_mean,
                                            const arma::cube prior_beta_cov, const arma::mat& beta, const arma::cube& Wtildei, const arma::mat& Utildei,
                                            const arma::mat& iSigmae){

  // # This program updated beta2, and requires no Metropolis step. It also uses
  // # the fact that the design matrix for the food is the same for every
  // # repeat. This code would have to be changed is this were not true.
  // #
  // # There is another special feature about this problem, namely that since
  // # there is no energy, Sigmae is diagonal. This makes the step for beta2 a
  // # lot different.
  // #
  // # Input
  // # Xtildei        = design matrix
  // # mmi            = # of recalls (assumed same for all)
  // # cov_prior_beta = the prior covariance matrix for beta1
  // # beta           = the current values of all the betas
  // # Wtildei        = the latent and non-latent responses
  // # Utildei        = the latent person-specific effects
  // # iSigmae        = the 2x2 inverse of Sigmae. iSigmae(1,1) = 1

  arma::mat xx         = Xtildei.slice(0);
  arma::mat cc2        = inv(inv(prior_beta_cov.slice(2)) + (mmi * iSigmae(2,2) * (xx.t() * xx)));
  double mmbeta        = beta.n_rows;
  arma::mat cc1        = arma::zeros(mmbeta,1);
  cc1        = cc1 + (inv(prior_beta_cov.slice(2)) * prior_beta_mean.col(2));
  for( int jji=0; jji< mmi; jji ++){
    arma::mat w = Wtildei.slice(jji);
    cc1    = cc1 + (iSigmae(2,2)*(xx.t() * (w.col(2) - Utildei.col(2))));
    cc1    = cc1 + (iSigmae(0,2)*(xx.t() * (w.col(0) - (xx * beta.col(0)) - Utildei.col(0))));
    cc1    = cc1 + (iSigmae(1,2)*(xx.t() * (w.col(1) - (xx * beta.col(1)) - Utildei.col(1))));
  }
  //cc2 = (cc2 + cc2.t())/2;
  arma::mat t = arma::sqrtmat_sympd(cc2);
  arma::colvec beta3      = (cc2 * cc1) + (t * arma::randn(mmbeta,1));

  return(beta3);

}

// [[Rcpp::export]]
arma::mat ginverse_c(const arma::mat& z,const double lambda){

  //   Computes the inverse of the Box-Cox transformation, and makes sure that
  // the argument does not get negative
  //


  arma::mat x(z.n_rows, z.n_cols);
  arma::mat one = arma::ones(z.n_rows, z.n_cols);
  if (lambda == 0){
    x = arma::exp(z);
  }
  if (lambda > 0){
    arma::mat w = one + (lambda * z);

    for(int i=0; i < w.n_rows; i++){
      for(int j =0; j< w.n_cols; j++){
        if(w(i,j) < 0){
          w(i,j) = 0;
        }

        x(i,j) = pow(w(i,j) , pow(lambda, -1));
      }
    }


  }

  return(x);

}




// [[Rcpp::export]]
arma::mat backtransform_c(const double lambda,arma::mat& Xtildei,arma::vec& beta,
                          double sigmae,double mumu,double sigsig,arma::vec& Utildei,double& n){
  // # Compute the 9 point backtransformation for any component
  // #
  // # INPUT
  // # lambda      = the tranformation parameter
  // # Xtildei     = the design matrix
  // # beta        = the posterior mean
  // # sigmae      = standard deviation of epsilon
  // # mumu        = the transformed mean
  // # sigsig      = the transformed standard deviation
  // # Utildei     = the realized U
  // # n           = the sample size
  //
  // # set the abscissas and weights
  arma::vec x = {-2.1,-1.3,-0.8,-0.5, 0.00, 0.5, 0.8, 1.3, 2.1};
  arma::vec w = {0.063345, 0.080255, 0.070458, 0.159698, 0.252489,0.159698, 0.070458, 0.080255, 0.063345};

  int t = x.n_elem;
  arma::vec first_term = (Xtildei * beta + Utildei);
  arma::mat temp = arma::repmat(first_term, 1, t )
    +  arma::repmat(x.t(), n, 1) * sqrt(2) * sigmae; //get the terms in the innermost parentheses
  temp =  (mumu + sigsig * temp / sqrt(2));  // de-standardize
  temp = ginverse_c(temp,lambda);
  temp = arma::sum(repmat(w.t(), n,1)%temp, 1); //get n*1 vector,  eq A.5


  return(temp);

}



// 
// 
// // [[Rcpp::export]]
// Rcpp::List perform_MCMC_c(const double nMCMC,const double nthin,const double n, const double mmi,const arma::cube& Xtildei,
//   arma::mat& Utildei,arma::mat& beta,arma::mat& alpha,arma::mat& GGalpha, arma::vec didconsume,
//   arma::mat& prior_alpha_cov, arma::mat& prior_alpha_mean, arma::cube Wtildei,arma::mat& iSigmae,arma::mat& iSigmau, 
//   arma::vec& Wistar, double r, double theta,double s22,double s33,arma::mat& Sigmau,arma::mat& Sigmae, arma::mat& prior_Sigmau_mean,
//   double prior_Sigmau_doff, double prior_Sigmae_doff,arma::mat& prior_beta_mean,arma::cube& prior_beta_cov,
//   double update_beta1_var_ind,double lambda_rec_food, double lambda_rec_energy,double mumu,double sigsig,double mu_e,double sig_e, 
//   double mdesign,double rw_ind,arma::vec& beta1_accept_count, double a0_food,double a0_energy,double ndist,double ndim){
//                                    
//   std::cout<< "Start the MCMC" << std::endl;
// // ###########################################################################
// // # Initialize the MCMC traces
// // ###########################################################################
// arma::mat  r_trace(nMCMC,1), theta_trace(nMCMC,1), s22_trace(nMCMC,1), s33_trace(nMCMC,1), alpha_trace(nMCMC,GGalpha.n_cols);
// arma::mat never_trace(nMCMC,1), usual_intake_food_trace(n, ndist), usual_intake_energy_trace(n, ndist);
//   
// std::cout<<"Check:1"<<std::endl;
// usual_intake_food_trace.fill(arma::datum::nan);  
// usual_intake_energy_trace.fill(arma::datum::nan);
// 
// arma::cube Sigmae_trace(3,3,nMCMC), Sigmau_trace(3,3,nMCMC), beta_trace(mdesign,3,nMCMC);
// 
// arma::mat temp_mat;
// 
// std::cout<<"Check:2"<<std::endl;
//   
// for(double jjMCMC = 0; jjMCMC< nMCMC; jjMCMC ++){
//     if(mod(jjMCMC,500) == 0.0){
//       std::cout<< "iteration <- " << jjMCMC<< std::endl;
//     }
//     
// // ###########################################################################
// // # Update Ni. You create this for everyone.
// // ###########################################################################
// 
// std::cout<<"Check:"<<std::endl;
// 
//  Rcpp::List  update_Ni = update_Ni_with_covariates_c(Xtildei,beta,Utildei,alpha,
//                                             GGalpha,n,mmi,didconsume);
//  arma::vec     Ni = update_Ni["Ni"];
//  arma::vec     ppi = update_Ni["ppi"];
//  arma::vec    isnever(Ni.n_elem);  // Indicator of a never-consumer
//  for(int i = 0; i < Ni.n_elem; Ni ++){
//    if(Ni(i)< 0){
//      isnever(i) = 1;
//    }else{
//      isnever(i) = 0;
//    }
//  }
// // ###########################################################################
// // # Update alpha. In the following, the complete con, ditional for alpha is
// // # that is a truncated normal from the left at alpha_min, but with mean (cc2
// // # * cc1) and variance cc2.
// // ###########################################################################
//     arma::mat   xx     = Xtildei.slice(1);
//     double mmnn        =  xx.n_cols;
//     arma::mat  cc1     = (inv(prior_alpha_cov)*prior_alpha_mean) + GGalpha.t()*Ni;
//     arma::mat  cc2     = inv(GGalpha.t()*GGalpha + inv(prior_alpha_cov));
//     arma::mat  mujj    = cc2 % cc1;
//    // cc2 = (cc2 + cc2.t())/2;
//     arma::mat sijj     = arma::sqrtmat_sympd(cc2);
//     arma::mat  alpha   = mujj + sijj*arma::randn(GGalpha.n_cols);
// // ###########################################################################
// // # Update W1 and W2
// // ###########################################################################
//     double  numgen     = 5;
//     arma::cube Wtildeinew = gen_Wtildei_1foodplusenergy_never_c(Wtildei,beta,Xtildei,Utildei,n,
//                                                       iSigmae,Wistar,mmi,numgen);
//     arma::cube Wtildei    = Wtildeinew;
// 
// // ###########################################################################
// // # Calculate W-XB-U
// // ###########################################################################
//     arma::mat  tt(n,ndim);
//     tt.zeros();
//         for(int jj = 0; jj < ndim; jj++){
//           tt.col(jj) = (Xtildei.slice(jj) * beta.col(jj)) + Utildei.col(jj);
//         }
// 
//      arma::cube  y(n, ndim, mmi);
//         for(int i = 0; i < mmi;i++){
//             y.slice(i) = tt;
//           }
//      arma::cube     qq = Wtildei - y;
// // ###########################################################################
// // # Update iSigmae
// // ###########################################################################
//              double  rnew    = updated_parameter_r_never_c(r,theta,s22,s33,qq,mmi,n);
//              double r        = rnew;
//              double thetanew = updated_parameter_theta_never_c(r,theta,s22,s33,qq,mmi);
//              double theta    = thetanew;
//              double s22new   = updated_parameter_s22_never_c(r,theta,s22,s33,qq,mmi,n);
//              double s22      = s22new;
//              double s33new   = updated_parameter_s33_never_c(r,theta,s22,s33,qq,mmi,n);
//              double s33      = s33new;
// 
//              arma::mat R = {{1, 0, r*cos(theta)}, 
//              {0,  1, r*sin(theta)},
//              {r*cos(theta), r*sin(theta), 1}};
//              arma::vec v = {1, sqrt(s22), sqrt(s33)};
//              arma::mat A = arma::diagmat(v);
//              
//              arma::mat Sigmae  = A * R * A; 
//              arma::mat iSigmae = inv(Sigmae);
// // ###########################################################################
// // # Update iSigmaU
// // ###########################################################################
//             Rcpp::List update_sig = update_iSigmau_c(Sigmau, prior_Sigmau_doff,
//                                            prior_Sigmau_mean,Utildei,n, jjMCMC);
// 
//             arma::mat Sigmau_new  = update_sig["Sigmau_new"];
//             arma::mat iSigmau_new = update_sig["iSigmau_new"];
//             arma::mat Sigmau      = Sigmau_new;
//               arma::mat iSigmau     = iSigmau_new;
// // ###########################################################################
// // # Update Utildei. This is done in two steps. In the first step, we generate
// // # it assuming that everyone is a consumer. In the second step, those who
// // # are never consumers, i.e., Ni < 0, have their values updated by a
// // # Metropolis step.
// // ###########################################################################
//             arma::mat Utildei_new= update_Utildei_c(Utildei,beta,Wtildei,iSigmae,
//                                            isnever,didconsume,Xtildei,mmi,iSigmau,n);
//             Utildei     = Utildei_new;
// //###########################################################################
// // # Update beta1 using a Metropolis Step.
// // ###########################################################################
//             double p = beta.n_rows;
//             arma::vec beta1(p), beta2(p), beta3(p);
//             if (rw_ind == 1){
//               beta1 = update_beta1_with_prior_mean_random_walk_c(Xtildei,mmi,
//                                                                   prior_beta_mean, prior_beta_cov,beta,Wtildei, Utildei,
//                                                                   iSigmae,isnever,update_beta1_var_ind);
//             } else {
//               beta1  = update_beta1_with_prior_mean_c(Xtildei,mmi,prior_beta_mean,
//                                                       prior_beta_cov,beta,Wtildei, Utildei,iSigmae,isnever,
//                                                       update_beta1_var_ind);
//             }
// // count if beta1 moves
//             arma::uvec accept_count = arma::find(beta1 == beta.col(1));
//             arma::vec beta1_accept_count = beta1_accept_count + (1 - accept_count);
//             beta.col(1)  = beta1;
// // ###########################################################################
// // # Update beta2. This does not need a Metropolis step
// // ###########################################################################
//               beta2 = update_beta2_with_prior_mean_c(Xtildei,mmi, prior_beta_mean,
//                                                       prior_beta_cov,beta,Wtildei, Utildei,iSigmae);
//               beta.col(2)  = beta2;
// // ###########################################################################
// // # Update beta2. This does not need a Metropolis step
// // ###########################################################################
//               beta3 = update_beta3_with_prior_mean_c(Xtildei,mmi, prior_beta_mean,
//                                                       prior_beta_cov,beta,Wtildei, Utildei,iSigmae);
//               beta.col(3)  = beta3;
// //###########################################################################
// // # Store results
// // ###########################################################################
//               Sigmae_trace.slice(jjMCMC) = Sigmae;
//               Sigmau_trace.slice(jjMCMC) = Sigmau;
//               beta_trace.slice(jjMCMC)   = beta;
//               r_trace(jjMCMC,1)          = r;
//               theta_trace(jjMCMC,1)      = theta;
//               s22_trace(jjMCMC,1)        = s22;
//               s33_trace(jjMCMC,1)        = s33;
//               alpha_trace.row(jjMCMC)      = alpha;
//               never_trace(jjMCMC,1)      = 1 - arma::accu(arma::normcdf(GGalpha * alpha))/n ;
// // ###########################################################################
// // # Compute distribution of usual intake.
// // # Suppose we have finished an MCMC step. In this step, we know who are
// // # non-consumers (N_i < 0), and who are consumers (N_i > 0).
// // # Use Gauss-Hermite quadrature method to approximate the Q_F,
// // # which is average amount of food on consumption day for consumers
// // # (equations A.5 in section A.16). This is done using the
// // # backtransform_20130925 function.
// // # Then plug it in to compute the usual intake for consumers (equation A.2
// // # in section A.15).
// // # Do this for about 200 MCMC steps near the end, with thinning of 50.
// // ###########################################################################
//               // double vt = (nMCMC - (ndist - 1) * nthin);
//               // arma::vec ut = arma::regspace(vt, nthin, nMCMC);
//               // 
//               // if(any( ut == jjMCMC)){
//               //   
//               //   arma::vec    uuindex(Ni.n_elem);  // Indicator of a never-consumer
//               //   for(int i = 0; i < Ni.n_elem; Ni ++){
//               //     if(Ni(i)< 0){
//               //       uuindex(i) = 1.0;
//               //     }else{
//               //       uuindex(i) = 0.0;
//               //     }
//               //   }
//               //   
//               //   arma::mat f = arma::sqrtmat_sympd(Sigmau);
//               //   arma::mat Utildei1 = arma::randn(n,Sigmau.n_cols)*f; 
//               //   arma::uvec row =  arma::find(uuindex == 1) ;   
//               //   double nindex = row.n_rows;
//               //   arma::mat temp_1[row.n_rows];
//               //   
//               //   arma::mat temp_mat = Xtildei.slice(2);
//               //   arma::vec temp_vec = Utildei.col(2);
//               //   temp_vec = temp_vec.elem(row);
//               //   temp_1 = backtransform_c(lambda_rec_food, temp_mat.rows(row) , beta.col(2), std::sqrt(Sigmae(2,2)),
//               //                          mumu, 
//               //                          sigsig, 
//               //                          temp_vec,
//               //                          row.n_rows);
//               //   // temp[temp < a0_food] = a0_food
//                 // temp = temp * pnorm(Xtildei[uuindex, ,1] %*% beta[ ,1]
//                 //                        + Utildei[uuindex,1])# get usual intake (n*1), eq A.2
//                 // usual_intake_food = temp
//                 // temp = backtransform_c(as.numeric(lambda_rec_energy), Xtildei[uuindex == 1, ,3],
//                 //                         beta[ ,3], sqrt(Sigmae[3,3]), mu_e, sig_e, Utildei[uuindex == 1,3],
//                 //                                                                           nindex)
//                 // temp[temp < a0_energy] = a0_energy
//                 // usual_intake_energy = temp
// // # store the results for this run
// //                 col_index = as.numeric((jjMCMC - (nMCMC - (ndist - 1)%*% nthin))/ nthin + 1)
// //                 usual_intake_food_trace[ ,col_index]= c(usual_intake_food, rep(NaN,(n - nindex)))
// //                 usual_intake_energy_trace[ ,col_index]=c(usual_intake_energy, rep(NaN,(n - nindex)))
// //                 
// //                 
// //               }
// //   }
// //   
// // ###########################################################################
// // # end of MCMC
// // ###########################################################################
//   std::cout<<'MCMC completed'<<std::endl;
//     std::cout<<' '<<std::endl;
//     std::cout<<' '<<std::endl;
//     std::cout<<' '<<std::endl;
//     std::cout<<' '<<std::endl;
//     std::cout<<' '<<std::endl;
//     std::cout<<' '<<std::endl;
// 
//    return Rcpp::List::create( Rcpp::Named("alpha_trace") = alpha_trace, Rcpp::Named("beta_trace") = beta_trace, Rcpp::Named("never_trace") = never_trace,
//                         Rcpp::Named("r_trace") = r_trace, Rcpp::Named("theta_trace") = theta_trace, Rcpp::Named("s22_trace") = s22_trace,
//                         Rcpp::Named("s33_trace") = s33_trace,Rcpp::Named("Sigmae_trace") = Sigmae_trace, Rcpp::Named("Sigmau_trace") =Sigmau_trace);  /*, usual_intake_food_trace = usual_intake_food_trace, usual_intake_energy_trace = usual_intake_energy_trace))*/
// //return 0;
// 
// }
// }