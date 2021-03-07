update_iSigmau <- function(Sigmau, prior_Sigmau_doff, prior_Sigmau_mean,Utildei,n, jjMCMC){
#library(MCMCpack)
    
aa         <- (prior_Sigmau_doff- nrow(Sigmau)-1)*(prior_Sigmau_mean)+crossprod(Utildei)
bb         <- prior_Sigmau_doff + n
#aa  = aa + (.0011 .* (1 .* n) .* eye(3));#ab 3/11
#bb  = bb + (.0011 .* (1 .* n));#ab 3/11
aa  <- aa + (.00001 * (1 * n) * pracma::eye(3))
bb  <- bb + (.00001 * (1 * n))
aa  <- (aa + t(aa)) / 2 #'

if (min(eigen(aa / bb)$values) <= 0.001){
    thecount2 <- thecount2 + 1
    if (thecount2 == 1){
    print(paste('Problem with Sigmau at step ',jjMCMC))
    print('Sigmau')
    print(Sigmau)
    print('Sigmae')
    print(Sigmae)
    print('aa / bb')
    aa / bb
    }
}
if (min(eigen(aa / bb)$values) > 0.001){
    Sigmau_new     <- MCMCpack::riwish(bb, aa)
    iSigmau_new    <- MASS::ginv(Sigmau_new)
}

return(list(Sigmau_new = Sigmau_new, iSigmau_new = iSigmau_new))

}
