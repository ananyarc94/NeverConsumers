update_Ni_with_covariates <- function(Xtildei,beta,Utildei,alpha,GGalpha,n,mmi,
                                      didconsume){
#
# This program updates the random variables Ni
#
# INPUT
# Xtildei        = design matrix
# mmi            = # of recalls (assumed same for all)
# beta           = the current values of all the betas
# Wtildei        = the latent and non-latent responses
# Utildei        = the latent person-specific effects
# alpha          = alpha in the never-consumer model
# didconsume     = Indicator that the person consumed on one of the
#                     recalls
tt <- (Xtildei[ , ,1]%*% beta[ ,1]) + Utildei[ ,1]
aq1 <- pnorm(tt,0,1)
aq2 <- pnorm(GGalpha%*%alpha,0,1)
cc1 <- aq2 * ( (1 - aq1) ^ mmi) / (1 - aq2)#ab 3/5/11
ppi <- 1 / (1 + cc1)# This is p_i in the paper
Ni      <- matrix(0,n,1)
#genww1  = gen_truncated_normals(-alpha,-alpha .* ones(n,1),50);
#genww2  = gen_truncated_normals(alpha,-alpha .* ones(n,1),50);
kk <- GGalpha %*% alpha#ab 3/11
genww1  <- gen_truncated_normals_never(-kk,-kk * pracma::ones(n,1),50)#ab 3/11
genww2  <- gen_truncated_normals_never(kk,-kk * pracma::ones(n,1),50)#ab 3/11

rri     <- as.numeric(pracma::rand(n,1) < ppi)
Ni      <- GGalpha %*% alpha + (didconsume * genww1) +
          ((1-didconsume) * (((1-rri) * genww1) - (rri * genww2)))
return(list(Ni = Ni, ppi = ppi))
}
