update_beta2_with_prior_mean <- function(Xtildei,mmi, prior_beta_mean,
                                         prior_beta_cov,beta,Wtildei, Utildei,iSigmae){

# This program updated beta2, and requires no Metropolis step. It also uses
# the fact that the design matrix for the food is the same for every
# repeat. This code would have to be changed is this were not true.
#
# There is another special feature about this problem, namely that since
# there is no energy, Sigmae is diagonal. This makes the step for beta2 a
# lot different.
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
cc2        <- MASS::ginv(MASS::ginv(prior_beta_cov[ , ,2]) +
                   (mmi * iSigmae[2,2] * (t(xx) %*% xx)))#
mmbeta     <- pracma::size(beta,1)
cc1        <- numeric(mmbeta)
cc1        <- cc1 + (MASS::ginv(prior_beta_cov[ , ,2]) %*% prior_beta_mean[ ,2])
for (jji    in 1:mmi){
    cc1    <- cc1 + (iSigmae[2,2] * (t(xx) %*% (Wtildei[ ,2,jji] - Utildei[ ,2])))
    cc1    <- cc1 + (iSigmae[1,2]* (t(xx) %*%
                    (Wtildei[ ,1,jji] - (xx %*% beta[ ,1]) - Utildei[ ,1])))
    cc1    <- cc1 + (iSigmae[2,3]* (t(xx) %*%
                    (Wtildei[ ,3,jji] - (xx %*% beta[ ,3]) - Utildei[ ,3])))
}
beta2      <- (cc2 %*% cc1) + (pracma::sqrtm(cc2)$B %*% pracma::randn(mmbeta,1))
return(beta2)
}
