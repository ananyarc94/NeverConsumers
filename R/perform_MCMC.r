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
    update_Ni <- update_Ni_with_covariates_c(Xtildei,beta,Utildei,alpha,
                                            GGalpha,n,mmi,didconsume)
    Ni       <- update_Ni$Ni
    ppi      <- update_Ni$ppi
    isnever  <- as.numeric(Ni < 0)# Indicator of a never-consumer
    
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
   
    ###########################################################################
    # Update W1 and W2
    ###########################################################################
    numgen     <- 5
    Wtildeinew <- gen_Wtildei_1foodplusenergy_never_c(Wtildei,beta,Xtildei,Utildei,n,
                                                    iSigmae,Wistar,mmi,numgen)
    Wtildei    <- round(Wtildeinew,4)
    
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
    #print(paste("iteration", i))
    rnew   <- updated_parameter_r_never_c(r,theta,s22,s33,qq,mmi,n)
    (r        <- rnew)
    #print(paste("iteration", i))
    #print(paste("r = ",r))
    thetanew <- updated_parameter_theta_never_c(r,theta,s22,s33,qq,mmi)
    (theta    <- thetanew)
    #print(paste("theta = ",theta))
    s22new   <- updated_parameter_s22_never_c(r,theta,s22,s33,qq,mmi,n)
    (s22      <- s22new)
    #print(paste("s22 = ",s22))
    s33new   <- updated_parameter_s33_never_c(r,theta,s22,s33,qq,mmi,n)
    (s33      <- s33new)
    #print(paste("s33 = ",33))

    R     <-  matrix(c(1,     0,    r*cos(theta),
                       0,     1,    r*sin(theta),
                       r*cos(theta), r*sin(theta), 1 ), nrow = 3, byrow = T)

    A <- diag(c(1, sqrt(s22), sqrt(s33)))
    Sigmae       <- A %*% R %*% A
    #print(paste("Sigmae = ", Sigmae))
    
    iSigmae      <- round(pracma::inv(Sigmae),3)
    #print(iSigmae)
    ###########################################################################
    # Update iSigmaU
    ###########################################################################
    
    update_sig <- update_iSigmau_c(Sigmau, prior_Sigmau_doff,
                                   prior_Sigmau_mean,Utildei,n,jjMCMC)

    Sigmau_new  <- update_sig$Sigmau_new
    iSigmau_new <- update_sig$iSigmau_new
    Sigmau      <- Sigmau_new

   # # # print(Sigmau)
    iSigmau     <- iSigmau_new
    ###########################################################################
    # Update Utildei. This is done in two steps. In the first step, we generate
    # it assuming that everyone is a consumer. In the second step, those who
    # are never consumers, i.e., Ni < 0, have their values updated by a
    # Metropolis step.
    ###########################################################################
    #print(iSigmae)
    #print(iSigmau)
    Utildei_new <- update_Utildei_c(Utildei,beta,Wtildei,iSigmae,
                                 Ni,isnever,didconsume,Xtildei,mmi,iSigmau,n)
    Utildei     <- Utildei_new
    
    
    # ###########################################################################
    # Update beta1 using a Metropolis Step.
    ###########################################################################
      if (rw_ind == 1) {
      beta1 <- update_beta1_with_prior_mean_random_walk_c(Xtildei,mmi,
                                                          prior_beta_mean, prior_beta_cov,beta,Wtildei, Utildei,
                                                          iSigmae,isnever,update_beta1_var_ind)
    } else {
      beta1 <- update_beta1_with_prior_mean_c(Xtildei,mmi,prior_beta_mean,
                                              prior_beta_cov,beta,Wtildei, Utildei,iSigmae,isnever,
                                              update_beta1_var_ind)
    }
    # count if beta1 moves
    beta1_accept_count <- beta1_accept_count + (1 - all(beta1 == beta[ ,1]))
    beta[ ,1]  <- beta1

    ###########################################################################
    # Update beta2. This does not need a Metropolis step
    ###########################################################################
    beta2 <- update_beta2_with_prior_mean_c(Xtildei,mmi, prior_beta_mean,
                                            prior_beta_cov,beta,Wtildei, Utildei,iSigmae)
    beta[ , 2]  <- beta2
   
    ###########################################################################
    # Update beta2. This does not need a Metropolis step
    ###########################################################################
    beta3 <- update_beta3_with_prior_mean_c(Xtildei,mmi, prior_beta_mean,
                                            prior_beta_cov,beta,Wtildei, Utildei,iSigmae)
    beta[ ,3]  <- beta3
    
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
    temp <- backtransform_c(as.numeric(lambda_rec_food), Xtildei[uuindex == 1, ,2],
                               beta[ ,2], sqrt(Sigmae[2,2]), mumu, sigsig,
                               Utildei1[uuindex == 1,2], nindex)
    temp[temp < a0_food] <- a0_food
    #print("Done upto usual_intake_food")
    temp <- temp * pnorm(Xtildei[uuindex == 1, ,1] %*% beta[ ,1]
                            + Utildei1[uuindex == 1,1])# get usual intake (n*1), eq A.2
    usual_intake_food <- as.vector(temp)
    temp <- backtransform_c(as.numeric(lambda_rec_energy), Xtildei[uuindex == 1, ,3],
                               beta[ ,3], sqrt(Sigmae[3,3]), mu_e, sig_e, Utildei1[uuindex == 1,3],
                               nindex)
    temp[temp < a0_energy] <- a0_energy
    usual_intake_energy <- as.vector(temp)
    #print("Done upto usual_intake_energy")
    # store the results for this run
    col_index = as.numeric((jjMCMC - (nMCMC - (ndist - 1) * nthin)) / nthin + 1)
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
