generate_usual_intake <- function(Sigmau_postmean, Sigmae_postmean,
                    beta_postmean, alpha_postmean, Xtildei, GGalpha, mumu,
                    sigsig,mu_e,sig_e, a0_food, a0_energy, lambda_rec_food,
                    lambda_rec_energy, B){
    ###########################################################################
# This program generate B realizations of usual intake.
###########################################################################
# INPUT:
#	Sigmau:             ture or posterior estimate of Sigmau
#	Sigmae:             true or posterior estimate of Sigmae
#	beta:               ture or posterior estimate of beta
#   alpha:              ture or posterior estimate of alpha
#   Xtildei:            covariates in the consumer model (n by p_beta by 3)
#   GGalpha:            covariates in the never-consumer model (n by p_alpha)
#   mumu:               mean of the food
#   sigsig:             standard deviation of the food
#   mu_e:               mean of energy
#   sig_e:              standard deviation of energy
#   a0_food:            half of the minimum of non-negative food, needed in back-transformation
#   a0_energy           half of the minimum of energy, needed in back-transformation
#   lambda_rec_food:    box-cox transformation parameter of food
#   lambda_rec_energy:  box-cox transformation parameter of energy
#   B:                  the number of 24HR recalls you want to generate
#
# OUTPUT:
#   data_wide:          generated data set with columns of the following order:
#                       recall food 1, ..., recall food B, recall energy 1, ... , recall energy B
###########################################################################
# Initialize the matrix for B realizations of usual intake.
###########################################################################
#   n: the sample size you want to generate, which is the same as the
#   number of rows of Xtildei.
 n <- pracma::size(Xtildei, 1)
 W <- array(NaN, c(n, 3, B))

###########################################################################
# Generate person-specific variation (only once)
###########################################################################
Utildei <- pracma::randn(n,size(Sigmau,2)) %*% pracma::sqrtm(Sigmau)$B
###########################################################################
# Generate the latent variable to model the probability of being consumer
###########################################################################
Ni      <- GGalpha%*%(alpha) + pracma::randn(n,1)

###########################################################################
# This is fixed part (not include day-to-day variation) of the 3-part model (equation 2).
###########################################################################
temp <- matrix(NaN,n, 3)
for (jj in 1:3){
    temp[ ,jj] <- (Xtildei[ , ,jj] %*% beta[ ,jj]) + Utildei[ ,jj]
}
###########################################################################
# Now get the B realizations of usual intake.
###########################################################################
for (b in 1:B){
     W[ , , b]  <- temp + (pracma::randn(n,3) %*% pracma::sqrtm(Sigmae)$B)# This is the right hand size of the 3-part model (equation 2)
}

 W[ ,1, ] <- (W[ ,1, ] > 0)
 W[ ,2, ] <-  (mumu + sigsig *  W[ ,2, ] / sqrt(2))
 # max(a0,(1 + (lambda .*     W(:,2,:)             ))) .^ (1 ./ lambda));
 W[ ,2, ] <- max(a0_food, ginverse(W[ ,2, ] , lambda_rec_food))

 W[ ,3, ] <-  (mu_e + sig_e *  W[ ,3, ] / sqrt(2))
 W[ ,3, ] <- max(a0_energy, ginverse(W[ ,3, ] , lambda_rec_energy))

 usual_intake_food <- matrix(as.numeric(Ni > 0), nrow = n, ncol = B, byrow = T) %*% t(W[ ,1, ]) %*% (W[ ,2, ])
 usual_intake_energy <- (W[ ,3, ])
 data_wide <- rbind(usual_intake_food, usual_intake_energy)
 return(data_wide)

}
