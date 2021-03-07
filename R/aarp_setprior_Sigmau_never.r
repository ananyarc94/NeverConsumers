aarp_setprior_Sigmau_never <- function(prior_Sigmau_mean){

Sigmau             <- prior_Sigmau_mean
for (jj in 1:3){
    for (kk in 1:3){
        if (abs(jj-kk) > 0){
            Sigmau[jj,kk] <- 0.50 * sqrt(Sigmau[jj,jj] * Sigmau[kk,kk])
        }
    }
}
prior_Sigmau_doff <- 1 + 1 + nrow(Sigmau)
prior_Sigmau_mean <- Sigmau

return(list(prior_Sigmau_mean = prior_Sigmau_mean, prior_Sigmau_doff = prior_Sigmau_doff))
}

