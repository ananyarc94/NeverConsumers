update_Utildei <- function(UtildeiNI,betaNI,WtildeiNI,iSigmaeNI,
                           NI_U,isneverNI,didconsumeNI,XtildeiNI,mmiNI,iSigmauNI,n){


###########################################################################
# Update UtildeiNI. This is done in two steps. In the first step, we generate
# it assuming that everyone is a consumer. In the second step, those who
# are never consumers, i.e., Ni < 0, have their values updated by a
# Metropolis step.
###########################################################################
#
# INPUT
# UtildeiNI        = the latent person-specific effects
# betaNI           = the current values of all the betaNIs
# WtildeiNI        = the latent and non-latent responses
# iSigmaeNI        = inverse of Sigmae
# NI_U             = the crucial random variable
# isneverNI        = Indicator that Ni < 0 and the person is a never consumer
# didconsumeNI     = Indicator that the person consumed on one of the
#                     recalls
# XtildeiNI        = design matrix
# mmiNI            = # of recalls (assumed same for all)
# iSigmauNI        = inverse of Sigmaeu


# Here is Step 1
UtildeiNI_Current <- UtildeiNI
qq <- (XtildeiNI[ , ,1] %*% betaNI[ ,1])
for (jj in 2:3){
    qq <- cbind(qq,  (XtildeiNI[ , ,jj] %*% betaNI[ ,jj]))
}


y = array(1, c(nrow(qq), ncol(qq), mmiNI))
for( i in 1:mmiNI){
  y[ , , i] = qq
}
ss <- WtildeiNI - y


c1  <- t(iSigmaeNI%*% t(rowSums(ss, dims = 2))) 
c2  <- pracma::inv(iSigmauNI + (mmiNI * iSigmaeNI))
#c2  = (c2+ c2')./2;
# Here is Step 2
isneverNI <- (NI_U < 0)

UtildeiNI <- t(c2 %*% t(c1)) + (pracma::randn(n,3)  %*% pracma::sqrtm(c2)$B)
c2_cand <- 1 / ((1 - pnorm(qq[ ,1] + UtildeiNI[ ,1])) ^ mmiNI)
c2_curr <- 1 / ( (1 - pnorm(qq[ ,1] + UtildeiNI_Current[ ,1])) ^ mmiNI)
rri     <- (pracma::rand(n,1) < (c2_cand / c2_curr))
ddff       <- pracma::size(UtildeiNI,2)
UtildeiNI_MH <- (UtildeiNI * (rri %*% pracma::ones(1,ddff)) ) +
             (UtildeiNI_Current * ((1 - rri) %*% pracma::ones(1,ddff)))

UtildeiNI_new <- (UtildeiNI * ((1 - isneverNI) %*% pracma::ones(1,ddff))) +
                (UtildeiNI_MH * (isneverNI %*% pracma::ones(1,ddff)))

for(i in 1:nrow(UtildeiNI_new)){
  for( j in 1:ncol(UtildeiNI_new)){
    UtildeiNI_new[i, j] = ifelse(is.na(UtildeiNI_new[i, j]), 0,UtildeiNI_new[i, j] )
  }
} 
return(UtildeiNI_new)
}

