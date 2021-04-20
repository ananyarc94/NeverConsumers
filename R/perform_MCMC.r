perform_MCMC = function(nMCMC,nthin, n, mmi, Xtildei, Utildei, beta, alpha, GGalpha, didconsume, 
                        prior_alpha_cov, prior_alpha_mean, Wtildei, iSigmae, iSigmau, Wistar,r, theta,
                        s22, s33, Sigmau, Sigmae, prior_Sigmau_mean, prior_Sigmau_doff, prior_Sigmae_doff,
                        prior_beta_mean, prior_beta_cov, update_beta1_var_ind, lambda_rec_food, 
                        lambda_rec_energy, mumu, sigsig, mu_e, sig_e, mdesign, rw_ind, beta1_accept_count, 
                        a0_food, a0_energy,ndist,ndim){
  
  print('Start the MCMC')
  ###########################################################################
  # Initialize the MCMC traces
  ###########################################################################
  r_trace      <- matrix(0,nMCMC,1)
  theta_trace  <- matrix(0,nMCMC,1)
  s22_trace    <- matrix(0,nMCMC,1)
  s33_trace    <- matrix(0,nMCMC,1)
  Sigmae_trace <- array(0,c(3,3,nMCMC))
  Sigmau_trace <- array(0,c(3,3,nMCMC))
  beta_trace   <- array(0,c(mdesign,3,nMCMC))
  alpha_trace  <- matrix(0,nMCMC,ncol(GGalpha))
  never_trace  <- matrix(0,nMCMC,1)   
  usual_intake_food_trace <- matrix(NaN, n, ndist)
  usual_intake_energy_trace <- matrix(NaN, n, ndist)
  
   for (jjMCMC in 1:nMCMC) {
      if (jjMCMC %% 500 == 0) {
        print(paste("iteration <- ", jjMCMC))
     }

   #print(paste("iteration <- ", jjMCMC))
    ###########################################################################
    # Update Ni. You create this for everyone.
    ###########################################################################
    update_Ni <- update_Ni_with_covariates(Xtildei,beta,Utildei,alpha,
                                            GGalpha,n,mmi,didconsume)
    Ni       <- update_Ni$Ni
    #sum(Ni)
    # Ni = c(1.2109,1.5063, -0.1293, -0.9696, 0.0762, 1.9260, 1.2277, 1.1190, 0.8241,
    #       1.8091, 0.8225, 2.7636, 0.8965, 1.9034, -0.0354, 0.3180, 0.3111, 0.4309,
    #      0.8641, 1.7307)
    ppi      <- update_Ni$ppi
    isnever  <- as.numeric(Ni < 0)# Indicator of a never-consumer
    
    # update_Ni<- update_Ni_with_covariates_c(Xtildei,beta,Utildei,alpha,
    #                                       GGalpha,n,mmi,didconsume)
    # Ni       <- update_Ni$Ni
    # ppi      <- update_Ni$ppi
    # isnever  <- as.numeric(Ni < 0)
    ###########################################################################
    # Update alpha. In the following, the complete con, ditional for alpha is
    # that is a truncated normal from the left at alpha_min, but with mean (cc2
    # * cc1) and variance cc2.
    ###########################################################################
    xx      <- Xtildei[ , ,1]
    mmnn    <- ncol(xx)
    cc1     <- (pracma::inv(prior_alpha_cov) %*% prior_alpha_mean) + crossprod(GGalpha, Ni)
    cc2     <- pracma::inv(t(GGalpha) %*% GGalpha + pracma::inv(prior_alpha_cov))
    mujj    <- cc2 %*% cc1
    sijj    <- pracma::sqrtm(cc2)$B
    alpha   <- mujj + sijj %*% pracma::randn(ncol(GGalpha),1)
    #alpha = c(1.3617,   -0.0815)
    ###########################################################################
    # Update W1 and W2
    ###########################################################################
    numgen     <- 5
    Wtildeinew <- gen_Wtildei_1foodplusenergy_never(Wtildei,beta,Xtildei,Utildei,n,
                                                    iSigmae,Wistar,mmi,numgen)
    Wtildei    <- round(Wtildeinew,4)
    #Wtildei = matrix(c(), nrow = 20, byrow = T)
    ###########################################################################
    # Calculate W-XB-U
    ###########################################################################
    tt <- matrix(0,n,ndim)
    for (jj in 1:ndim) {
      tt[ ,jj] <- (Xtildei[ , ,jj] %*% beta[ ,jj]) + Utildei[ ,jj]
    }
    
    y = array(1, c(nrow(tt), ncol(tt), mmi))
    for (i in 1:mmi) {
      y[ , , i] = tt
    }
    qq <- Wtildei - y
    ###########################################################################
    # Update iSigmae
    ###########################################################################
    rnew   <- updated_parameter_r_never(r,theta,s22,s33,qq,mmi,n)
    (r        <- rnew)
    thetanew <- updated_parameter_theta_never(r,theta,s22,s33,qq,mmi)
    (theta    <- thetanew)
    s22new   <- updated_parameter_s22_never(r,theta,s22,s33,qq,mmi,n)
    (s22      <- s22new)
    s33new   <- updated_parameter_s33_never(r,theta,s22,s33,qq,mmi,n)
    (s33      <- s33new)

    R     <-  matrix(c(1,     0,    r*cos(theta),
                       0,     1,    r*sin(theta),
                       r*cos(theta), r*sin(theta), 1 ), nrow = 3, byrow = T)

    A <- diag(c(1, sqrt(s22), sqrt(s33)))
    Sigmae       <- A %*% R %*% A
    #paste("Sigmae", Sigmae)
    
    #Sigmae = matrix(c(1.0000,         0,    0.0848,
    #                  0,    2.1192,    0.4257,
    #                  0.0848,    0.4257,    1.5144), nrow = 3, byrow = T)
    iSigmae      <- round(pracma::inv(Sigmae),3)
    #print(iSigmae)
    ###########################################################################
    # Update iSigmaU
    ###########################################################################
    
    update_sig <- update_iSigmau(Sigmau, prior_Sigmau_doff,
                                   prior_Sigmau_mean,Utildei,n,jjMCMC)

    Sigmau_new  <- update_sig$Sigmau_new
    iSigmau_new <- update_sig$iSigmau_new
    Sigmau      <- Sigmau_new

   #  Sigmau = Sigmau_new = matrix(c(0.4452,    0.2203,    0.0856,
   #                    0.2203,    0.5566,    0.0273,
   #                    0.0856,    0.0273,    0.5109), nrow = 3, byrow = T)
   # # # print(Sigmau)
    iSigmau     <- pracma::inv(Sigmau_new)
    ###########################################################################
    # Update Utildei. This is done in two steps. In the first step, we generate
    # it assuming that everyone is a consumer. In the second step, those who
    # are never consumers, i.e., Ni < 0, have their values updated by a
    # Metropolis step.
    ###########################################################################
    #print(iSigmae)
    #print(iSigmau)
    Utildei_new <- update_Utildei(Utildei,beta,Wtildei,iSigmae,
                                 Ni,isnever,didconsume,Xtildei,mmi,iSigmau,n)
    Utildei     <- Utildei_new
    
    # Utildei = matrix(c(0.7441,    0.6363,   -0.4027,
    #                    -0.9594,    0.0139,    0.6437,
    #                    0.0125,   -0.3568,    0.1962,
    #                    0.1495,    0.0328,    0.7694,
    #                    -0.4692,   -0.0193,   -0.3197,
    #                    -0.5966,   -0.3704,   -0.2199,
    #                    -0.1015,    0.5187,    0.1877,
    #                    -0.3132,   -0.4491,   -0.4452,
    #                    -0.4777 ,  -0.1950,   -0.1410,
    #                    -0.5459 ,  -0.0889 ,   0.7332,
    #                    0.9566 ,   0.3427 ,  -0.2191,
    #                    -0.4406 ,  -0.3301,   -0.8149,
    #                    0.9720  ,  0.8079 ,  -0.9089,
    #                    0.2645  , -0.0624 ,   0.5897,
    #                    -0.2668 ,   0.3419 ,   0.0790,
    #                    -0.0198 ,   0.0818 ,   0.0049,
    #                    0.4486  , -0.1085  , -0.1987,
    #                    1.2091  ,  0.5237  , -0.5848,
    #                    0.3839  ,  0.4808  , -0.4679,
    #                    -0.4724 ,  -0.7147 ,   0.8115), ncol = 3, byrow = T)
    # ###########################################################################
    # Update beta1 using a Metropolis Step.
    ###########################################################################
    
    # beta_post = matrix(c(-1.9560,   -0.4660 ,  -0.0020,
    #                      0.2860,    0.2774   ,-0.2089,
    #                      0.0075 ,   0.2914    ,0.3137,
    #                      1.6245  ,  0.5700    ,0.2958,
    #                      -0.6466  , -0.1751    ,0.2155), nrow = 5, byrow = T)
    # 
    # 
    
    if (rw_ind == 1) {
      beta1 <- update_beta1_with_prior_mean_random_walk(Xtildei,mmi,
                                                          prior_beta_mean, prior_beta_cov,beta,Wtildei, Utildei,
                                                          iSigmae,isnever,update_beta1_var_ind)
    } else {
      beta1 <- update_beta1_with_prior_mean(Xtildei,mmi,prior_beta_mean,
                                              prior_beta_cov,beta,Wtildei, Utildei,iSigmae,isnever,
                                              update_beta1_var_ind)
    }
    # count if beta1 moves
    beta1_accept_count <- beta1_accept_count + (1 - all(beta1 == beta[ ,1]))
    beta[ ,1]  <- beta1
    
    #beta[ ,1]  <- beta_post[ , 1]
    ###########################################################################
    # Update beta2. This does not need a Metropolis step
    ###########################################################################
    beta2 <- update_beta2_with_prior_mean(Xtildei,mmi, prior_beta_mean,
                                            prior_beta_cov,beta,Wtildei, Utildei,iSigmae)
    beta[ , 2]  <- beta2
    #beta[ ,2]  <- beta_post[ , 2]
    ###########################################################################
    # Update beta2. This does not need a Metropolis step
    ###########################################################################
    beta3 <- update_beta3_with_prior_mean(Xtildei,mmi, prior_beta_mean,
                                            prior_beta_cov,beta,Wtildei, Utildei,iSigmae)
    beta[ ,3]  <- beta3
    #beta[ ,3]  <- beta_post[ , 3]
    ###########################################################################
    # Store results
    ###########################################################################
    Sigmae_trace[ , ,jjMCMC] <- Sigmae
    Sigmau_trace[ , ,jjMCMC] <- Sigmau
    beta_trace[ , ,jjMCMC]   <- beta
    r_trace[jjMCMC,1]        <- r
    theta_trace[jjMCMC,1]    <- theta
    s22_trace[jjMCMC,1]      <- s22
    s33_trace[jjMCMC,1]      <- s33
    alpha_trace[jjMCMC, ]    <- alpha
    never_trace[jjMCMC,1]    <- 1 - sum(pnorm(GGalpha %*% alpha))/n
    
    
    ###########################################################################
    # Compute distribution of usual intake.
    # Suppose we have finished an MCMC step. In this step, we know who are
    # non-consumers (N_i < 0), and who are consumers (N_i > 0).
    # Use Gauss-Hermite quadrature method to approximate the Q_F,
    # which is average amount of food on consumption day for consumers
    # (equations A.5 in section A.16). This is done using the
    # backtransform_20130925 function.
    # Then plug it in to compute the usual intake for consumers (equation A.2
    # in section A.15).
    # Do this for about 200 MCMC steps near the end, with thinning of 50.
    ###########################################################################
    v = as.numeric(nMCMC - (ndist - 1) %*% nthin)
    ut = seq(v,nMCMC, by = nthin)
    if (any( ut == jjMCMC)) {
    uuindex <- as.numeric(Ni > 0)
    nindex <- sum(uuindex)
    Utildei1 <- pracma::randn(n,pracma::size(Sigmau,2)) %*% pracma::sqrtm(Sigmau)$B
    temp <- backtransform(as.numeric(lambda_rec_food), Xtildei[uuindex == 1, ,2],
                               beta[ ,2], sqrt(Sigmae[2,2]), mumu, sigsig,
                               Utildei1[uuindex == 1,2], nindex)
    temp[temp < a0_food] <- a0_food
    #print("Done upto usual_intake_food")
    temp <- temp * pnorm(Xtildei[uuindex, ,1] %*% beta[ ,1]
                            + Utildei1[uuindex,1])# get usual intake (n*1), eq A.2
    usual_intake_food <- as.vector(temp)
    temp <- backtransform(as.numeric(lambda_rec_energy), Xtildei[uuindex == 1, ,3],
                               beta[ ,3], sqrt(Sigmae[3,3]), mu_e, sig_e, Utildei1[uuindex == 1,3],
                               nindex)
    temp[temp < a0_energy] <- a0_energy
    usual_intake_energy <- as.vector(temp)
    #print("Done upto usual_intake_energy")
    # store the results for this run
    col_index = as.numeric((jjMCMC - (nMCMC - (ndist - 1) %*% nthin)) / nthin + 1)
    usual_intake_food_trace[ ,col_index] <- c(usual_intake_food, rep(NaN,(n - nindex)))
    usual_intake_energy_trace[ ,col_index] <- c(usual_intake_energy, rep(NaN,(n - nindex)))
    }
  }
  
  ##########################################################################
  # end of MCMC
  ###########################################################################
  print('MCMC completed')
  print(' ')
  print(' ')
  print(' ')
  print(' ')
  print(' ')
  print(' ')
  
  return(list(alpha_trace = alpha_trace, beta_trace = beta_trace, never_trace = never_trace,
              r_trace = r_trace, theta_trace = theta_trace, s22_trace = s22_trace, 
              s33_trace = s33_trace,Sigmae_trace = Sigmae_trace, Sigmau_trace = Sigmau_trace,
              usual_intake_food_trace = usual_intake_food_trace, usual_intake_energy_trace = usual_intake_energy_trace)) 
}
