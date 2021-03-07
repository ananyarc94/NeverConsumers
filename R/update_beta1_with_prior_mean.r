update_beta1_with_prior_mean <- function(Xtildei,mmi, prior_beta_mean,
                                    prior_beta_cov,beta,Wtildei, Utildei,
                                    iSigmae,isnever,update_beta1_var_ind){

# This program updated beta1 using  Metropolis step. It also uses the fact
# that the design matrix for the food is the same for every repeat. This
# code would have to be changed is this were not true.
#
# Input
# Xtildei        = design matrix
# mmi            = # of recalls (assumed same for all)
# cov_prior_beta = the prior covariance matrix for beta1
# beta           = the current values of all the betas
# Wtildei        = the latent and non-latent responses
# Utildei        = the latent person-specific effects
# iSigmae        = the 2x2 inverse of Sigmae. iSigmae(1,1) = 1
xx         <- (Xtildei[ , ,1])
cc2        <- MASS::ginv(MASS::ginv(prior_beta_cov[ , ,1]) +
                (mmi * iSigmae[1,1] * (crossprod(xx))))
mmbeta     <- nrow(beta)
cc1        <- matrix(0,mmbeta,1)
cc1        <- cc1 + (MASS::ginv(prior_beta_cov[ , ,1]) %*% prior_beta_mean[ ,1])
for (jji    in 1:mmi){
    cc1    <- cc1 + (iSigmae[1,1]* (t(xx) %*% (Wtildei[ ,1,jji] - Utildei[ ,1])))
    cc1    <- cc1 + (iSigmae[1,2] * (t(xx) %*%
               (Wtildei[ ,2,jji] - (xx %*% beta[ ,2]) - Utildei[ ,2])))
    cc1    <- cc1 + (iSigmae[1,3]* (t(xx) %*%
               (Wtildei[ ,3,jji] - (xx %*% beta[ ,3]) - Utildei[ ,3])))
}
beta1_cand <- (cc2 %*% cc1) + (pracma::sqrtm(cc2/update_beta1_var_ind)$B %*% pracma::randn(mmbeta,1))
beta1_curr <- beta[ ,1]
lc2_cand   <- -mmi * sum(is.finite(isnever *
               log((1 - pnorm((xx %*% beta1_cand) + Utildei[ ,1])))))
lc2_curr   <- -mmi * sum(is.finite(isnever *
               log((1 - pnorm((xx%*% beta1_curr) + Utildei[ ,1])))))
lc1_cand   <- t(cc1) %*% beta1_cand - t(beta1_cand) %*% MASS::ginv(cc2) %*% beta1_cand /2
lc1_curr   <- t(cc1) %*% beta1_curr - t(beta1_curr) %*% MASS::ginv(cc2) %*% beta1_curr /2

gghh       <- min(1,exp( (1-update_beta1_var_ind) %*% lc1_cand + lc2_cand
                      - (1-update_beta1_var_ind) %*% lc1_curr - lc2_curr))
rri        <- as.numeric((pracma::rand(1,1) < gghh))
beta1      <- (beta1_cand * rri) + (beta1_curr * (1 - rri))

return(beta1)

}
# $$$ beta1_curr = beta(:,1);
# $$$ beta1_cand = beta1_curr + (sqrtm(2.*cc2) * randn(mmbeta,1));
# $$$
# $$$ lc1_cand = cc1'*beta1_cand - 0.5 .* beta1_cand'*inv(cc2)*beta1_cand;
# $$$ lc1_curr = cc1'*beta1_curr - 0.5 .* beta1_curr'*inv(cc2)*beta1_curr;
# $$$ lc2_cand   = -mmi .* sum(isnever .* ...
# $$$                log((1 - normcdf((xx * beta1_cand) + Utildei(:,1)))));
# $$$ lc2_curr   = -mmi .* sum(isnever .* ...
# $$$                log((1 - normcdf((xx * beta1_curr) + Utildei(:, ...
# $$$                                                   1)))));
# $$$ gghh       = min(1,exp(lc2_cand + lc1_cand - lc2_curr -lc1_curr));
# $$$ rri        = (rand(1,1) < gghh);
# $$$ beta1      = (beta1_cand .* rri) + (beta1_curr .* (1 - rri));
