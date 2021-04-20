update_Utildei <- function(Utildei,beta,Wtildei,iSigmae,
                           Ni,isnever,didconsume,Xtildei,mmi,iSigmau,n){


###########################################################################
# Update Utildei. This is done in two steps. In the first step, we generate
# it assuming that everyone is a consumer. In the second step, those who
# are never consumers, i.e., Ni < 0, have their values updated by a
# Metropolis step.
###########################################################################
#
# INPUT
# Utildei        = the latent person-specific effects
# beta           = the current values of all the betas
# Wtildei        = the latent and non-latent responses
# iSigmae        = inverse of Sigmae
# Ni             = the crucial random variable
# isnever        = Indicator that Ni < 0 and the person is a never consumer
# didconsume     = Indicator that the person consumed on one of the
#                     recalls
# Xtildei        = design matrix
# mmi            = # of recalls (assumed same for all)
# iSigmau        = inverse of Sigmaeu


# Here is Step 1
Utildei_Current <- Utildei
qq <- (Xtildei[ , ,1] %*% beta[ ,1])
for (jj in 2:3){
    qq <- cbind(qq,  (Xtildei[ , ,jj] %*% beta[ ,jj]))
}


y = array(1, c(nrow(qq), ncol(qq), mmi))
for( i in 1:mmi){
  y[ , , i] = qq
}
ss <- Wtildei - y


c1  <- t(iSigmae%*% t(rowSums(ss, dims = 2))) 
c2  <- pracma::inv(iSigmau + (mmi * iSigmae))
#c2  = (c2+ c2')./2;
# Here is Step 2
isnever <- (Ni < 0)

Utildei <- t(c2 %*% t(c1)) + (pracma::randn(n,3)  %*% pracma::sqrtm(c2)$B)
c2_cand <- 1 / ((1 - pnorm(qq[ ,1] + Utildei[ ,1])) ^ mmi)
c2_curr <- 1 / ( (1 - pnorm(qq[ ,1] + Utildei_Current[ ,1])) ^ mmi)
rri     <- (pracma::rand(n,1) < (c2_cand / c2_curr))
ddff       <- pracma::size(Utildei,2)
Utildei_MH <- (Utildei * (rri %*% pracma::ones(1,ddff)) ) +
             (Utildei_Current * ((1 - rri) %*% pracma::ones(1,ddff)))

Utildei_new <- (Utildei * ((1 - isnever) %*% pracma::ones(1,ddff))) +
                (Utildei_MH * (isnever %*% pracma::ones(1,ddff)))

for(i in 1:nrow(Utildei_new)){
  for( j in 1:ncol(Utildei_new)){
    Utildei_new[i, j] = ifelse(is.na(Utildei_new[i, j]), 0,Utildei_new[i, j] )
  }
} 
return(Utildei_new)
}

