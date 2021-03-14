gen_Wtildei_1foodplusenergy_never <- function(Wtildei,beta,Xtildei,Utildei,n,
                                              iSigmae,Wistar,mmi,numgen){

# The method of construction is the following. In the notes, W_{ijk}
# denotes the ith person, the jth VARIABLE, and the kth REPLICATE.
#
# On the other hand, Xtilde(ii,:,jj) is the design matrix for the jth
# variable. Similarly, Utildei(:,:,jj) is the random effects for the jth
# variable. In this case, we are interested in the 2nd variable.

#library(iemisc)
Wtildeinew  <- Wtildei
Wtildeio    <- Wtildei
########################################################################
# Update W_{i1k}
########################################################################
varnum      <- 1# The variable you want to generate from the truncated normal
for (kk      in 1:mmi){# This is the replicate #ab 3/5/11
    C2      <- 1/ iSigmae[varnum,varnum]
    C1      <- iSigmae[varnum,varnum] * ((Xtildei[ , ,varnum] %*% beta[ ,varnum])
                                         + Utildei[ ,varnum])
    for (jj  in 1:pracma::size(Xtildei,3)){
        if (abs(jj - varnum) > 0){
            qq <- (Wtildei[ ,jj,kk] - (Xtildei[ , ,jj] %*% beta[ ,jj]) -
                       Utildei[ ,jj])
            # Note that this is W_{i,jj,kk} in the notes
            C1 <- C1 - (iSigmae[varnum,jj] * qq)
        }
    }
    mu      <- C2 * C1
    sigma   <- sqrt(C2)
    startxi <- mu/sigma
    genww1  <- gen_truncated_normals_never_c(-mu/sigma,-startxi,numgen)
    genww2  <- gen_truncated_normals_never_c(mu/sigma,-startxi,numgen)
    Wtildeinew[ ,varnum,kk] <- mu + (sigma * ((Wistar[ ,kk]* genww1)
                    - ((1 - Wistar[ ,kk]) * genww2)))
}
########################################################################
# Update W_{i2k}. You only do this for the cases that W_{i1k} < 0;
########################################################################
varnum      <- 2# The variable you want to generate from the truncated normal
for (kk      in 1:mmi){# This is the replicate #ab 3/5/11
    C2      <- 1 / iSigmae[varnum,varnum]
    C1      <- iSigmae[varnum,varnum] * ((Xtildei[ , ,varnum] %*% beta[ ,varnum])
                                          + Utildei[ ,varnum])
    for (jj  in 1:pracma::size(Xtildei,3)){
        if (abs(jj - varnum) > 0){
            qq <- (Wtildeinew[ ,jj,kk] - (Xtildei[ , ,jj] %*% beta[ ,jj])
                   - Utildei[ ,jj])
            # Note that this is W_{i,jj,kk} in the notes
            C1 <- C1 - (iSigmae[varnum,jj] * qq)
        }
    }
    mu      <- C2 * C1
    sigma   <- sqrt(C2)
    Wtildeinew[ ,varnum,kk] <- mu + (sigma * pracma::randn(n,1))
}

 # Wtildeinew[ ,varnum, ] <- ((Wtildei[ ,varnum, ]) * Wistar)
 #        + (Wtildeinew[ ,varnum, ] * (1-Wistar))
 # 
 a = ((Wtildei[ ,varnum, ]) * Wistar)
 b = ((Wtildeinew[ ,varnum, ]) * (1 - Wistar))

 Wtildeinew[ ,varnum, ] = a + b
 return(Wtildeinew = Wtildeinew)
}
