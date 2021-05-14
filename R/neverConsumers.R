#' Identify Never-Consumers in Zero Inflated Poisson Models
#'
#'@description A function to estimate percentage of never-consumers in the context of nutrition. The method implemented in this function has been developed in the paper
#'"Measurement error models with zero - inflation and multiple sources of zero with an application to the never consumers problem in nutrition". If we have data on whether 
#'or not a person has consumed a certain food on the day of the study, amount of food consumed and energy obtained with a bunch of covariates and multiple recalls then we 
#'can use this function to model the data and obtain an estimate of the percentage of never-consumers of the food of interest. The data for this function must include 
#'number of recalls and the response in three different variable (i.e. an indicator of food consumed, amount of food consumed and amount of energy generated from the food). 
#'It also needs some covariates. A list of other parameters that must be specified by the user are given below. 
#' 
#' @usage  neverConsumers(zz, mmi, X_cols, Food_col, Energy_col, ID_col, FFQ_food_col,
#' FFQ_energy_col, with_covariates_ind, n_gen, lambda_rec_food, lambda_rec_energy, 
#' lambda_FFQ_food, lambda_FFQ_energy,beta_temp, Sigmau_temp_episodically)
#' 
#' @param zz  The dataset in the form of a dataframe with headers.
#' @param mmi  Number of Recalls, it should be an integer.
#' @param X_cols  Column id's specifying where the covariates are
#' @param Food_col  Column id specifying where the episodic variable is
#' @param Energy_col  Column id specifying where energy (continuous) is
#' @param ID_col  Column id for ID of the subjects      
#' @param FFQ_food_col  Column where the FFQ for the food is. Set = [] if no FFQ. Please be aware that if there is a FFQ, 
#' you really should add this in, because FFQ have a tendency to generate massive, high leverage outliers, as in the EATS data
#' @param FFQ_energy_col  Column where the FFQ for energy is. Set = [] if no FFQ Please be aware that if there is a FFQ, 
#' you really should add this in, because FFQ have a tendency to generate massive, high leverage outliers, as in the EATS data
#' @param with_covariates_ind  What to include as covariates in the ever consumer model; 0: a column of ones; 1: a column of 
#' ones, the FFQ, and the indicator that the FFQ=0; 2: a column of ones and the FFQ; 3: a column of ones and the indicator 
#' that the FFQ=0. By default the value is set to 3.
#' @param n_gen  Number of realizations of usual intake to generate. Must be a positive integer, or 0 if no realizations to
#'  be generated. Default is 1.
#' @param lambda_rec_food  lambda value required for BoxCox transformation on food. Default is 0.5.
#' @param lambda_rec_energy  lambda value required for BoxCox transformation on energy. Default is 0.5.
#' @param lambda_FFQ_food  lambda value required for BoxCox transformation on FFQ for the food. Default is 0.5.
#' @param lambda_FFQ_energy  lambda value required for BoxCox transformation on FFQ for the energy. Default is 0.5.
#' @param beta_temp  initial value for beta_temp. By default it is a matrix of zero's.
#' @param Sigmau_temp_episodically  initial value for Sigmau_temp. By default a matrix with diagonal elements 1 and off-diagonals 0.5
#' @param nburn  Size of the burnin. Default value is 10000.
#' @param nMCMC  Number of MCMC iterations. The more the better. Default value is 50000.
#' @param nthin  Thining. Default value is 50.
#' @param ndist  average of the last ndist MCMC steps (after thinning) to get the cdf of usual intake. Default value is 200. 
#' @param beta_start_ind  the starting value of beta. By default set to 9999: episodically; all other number: a 5*3 matrix 
#' with the number as each cell. This should be an integer.
#' @param beta_prior_mean_ind  the prior mean for beta 9999: episodically (estimated by Saijuan's code) 0: regular 
#' (a 5*3 matrix with all elements = 0). Please specify either 9999 or 0. 
#' @param rw_ind  do you want to use the random walk proposal for beta_1? 1: yes, use random walk proposal with 
#' Normal(beta_1, C_2 * M);  0: no, use Normal(C_2*C_1, C_2*M). Please specify either 1 or 0. Default value is 1 
#' @param update_beta1_var_ind  the variance for updating beta1, i.e. this is the M in 
#' Normal(beta_1, C_2*M) and Normal(C_2 * C_1, C_2 *M). It has to be a positive number. Default value is 0.5
#' @param Sigmau_start_ind  the starting value of Sigmau; 1: episodically (estimated by Saijuan's code); 
#' 2: regular (a 3*3 matrix with diagonal elements = 1 and offdiagonal elements = 0.5) .Please specify either 1 or 2. Default value is 1
#' @param Sigmau_prior_mean_ind  the prior mean for Sigmau; 1: episodically (estimated by Saijuan's code); 
#' 2: regular (a 3*3 matrix with diagonal elements = 1 and off-diagonal elements = 0.5); 3: half of regular. Please specify either 1 or 2. 
#' Default value is 2
#' @param ndim  the number of dimensions, here by default set to 3 for indicator, amount and energy. Should be an integer.
#' @param beta1_accept_count  count how many times beta1 moves. Default value is 0.
#' @param myseed  initialize random seed
#'
#' @return 
#' alpha_postmean,alpha_postsd, alpha_ci      : Mean, Sd and CI's of the mean parameter of the prior distribution of percentage 
#'                                              of never consumers \cr
#' \cr                                              
#' never_postmean, never_postsd, never_ci     : Mean, Sd and CI's of the percentage of never consumers \cr 
#' \cr                      
#' beta_postmean, beta_postsd, beta_ci        : Mean, Sd and CI's of slope parameters for the covariates \cr
#' \cr                        
#' Sigmau_postmean, Sigmau_postsd, Sigmau_ci  : Mean, Sd and CI's of the variance covariance matrix of the random effects \cr   
#' \cr      
#' Sigmae_postmean, Sigmae_postsd, Sigmae_ci  : Mean, Sd and CI's of the variance covariance matrix of the white noise  \cr      
#' \cr      
#' mu_ui_food, sig_ui_food                    : Mean and Sd of the distribution of usual intake of food   \cr  
#' \cr                     
#' mu_ui_energy, sig_ui_energy                : Mean and Sd of the distribution of usual intake of the ratio = food/(energy/1000) \cr
#' \cr
#' mu_ui_ratio, sig_ui_ratio                  : Mean and Sd of the distribution of usual intake of energy \cr   
#' \cr                        
#' food_distribution, energy_distribution     : Distribution of usual intake of food, energy and food/(energy/1000)                        
#' and ratio_distribution
#' @export
#'           
#' @examples 
#' 
#' # load the data
#' zz <- read.csv(file.choos(), header = T)
#' 
#' # specify all the required parameters with no default values
#' X_cols          <- 2:5 # Columns where the covariates are
#' Food_col        <- 6 # Column where the episodic variable is
#' Energy_col      <- 7 # Column where energy (continuous) is
#' ID_col          <- 1 # Column for ID
#' FFQ_food_col    <- 4 # Column where the FFQ for the food is. 
#' FFQ_energy_col  <- 5 # Column where the FFQ for energy is. 
#' with_covariates_ind <- 3# What to include as covariates in the
#' mmi             <- 4# Number of recalls, integer
#' 
#' # make the function call as follows:
#' 
#' neverConsumers(zz, mmi, X_cols, Food_col, Energy_col, ID_col, FFQ_food_col,
#' FFQ_energy_col, with_covariates_ind)
#' 
#' 
neverConsumers = function(zz, mmi, X_cols, Food_col, Energy_col, ID_col, FFQ_food_col,
                          FFQ_energy_col, with_covariates_ind = 3, n_gen = 1, lambda_rec_food = 0.5, 
                          lambda_rec_energy = 0.5, lambda_FFQ_food = 0.5, lambda_FFQ_energy = 0.5, beta_temp = NULL,
                          Sigmau_temp_episodically = NULL, nburn = 10000, nMCMC = 50000, nthin = 50, ndist = 200, 
                          beta_start_ind = 9999, beta_prior_mean_ind = 0, rw_ind = 1, update_beta1_var_ind = 0.5, 
                          Sigmau_start_ind = 1, Sigmau_prior_mean_ind = 2, ndim = 3, 
                          beta1_accept_count = 0, myseed = 6309021, eps = 10^-4){
  #####################################################################################
  # Compatibility checks
  #####################################################################################
  if (ndist * nthin > nMCMC) {
    stop(paste('Please decrease ndist or increase nMCMC'))
  }
  if (!is.data.frame(zz)) {
    stop(paste("The data must be passed as a dataframe."))
  }
  if (is.null(mmi) || mmi <  1) {
    stop(paste("MMI (#recalls) needs to be specified and be > 1. It should also be an integer."))
  }
  if (mmi %% 1 != 0) {
    stop(paste("MMI should be an integer corresponding to the number of recalls in your dataset"))
  }
  if (nrow(zz) %% mmi != 0) {
    stop(paste("Number of recalls (MMI) doesn't match with the dataset"))
  }
  if (is.null(X_cols)) {
    stop(paste("Please specify where the column ID's for the covariates are"))
  }
  if (is.null(Food_col) || is.null(Energy_col)) {
    stop(paste("Please specify where the column ID's for energy and food values are"))
  }
  if (is.null(FFQ_food_col) || is.null(FFQ_energy_col)) {
    stop(paste("Please specify where the column ID's for food frequency questionnaire energy and food 
               values are. They can be among your covariates as well."))
  }
  if (is.null(ID_col)) {
    stop(paste("Please specify where the column with the subject ID's"))
  }
  # Initilization of SigmaU
  if (is.null(Sigmau_temp_episodically)) {
    Sigmau_temp_episodically = diag(rep(0.5,3)) + matrix(0.5, 3, 3)
  }
  
  ###########################################################################
  # Get the data. Assumes the first row are the variable names
  ###########################################################################
  
  zz         <-  zz[is.finite(rowSums(zz)), ]  ## get rid of rows with missing values
  ###########################################################################
  # How many subjects are there?
  ###########################################################################
  ID          = zz$ID
  n_subject   = length(unique(ID))
  ID = rep(1:n_subject, each = mmi)
  ###########################################################################
  # If you do not have an FFQ for the food, you need to set G. In this case,
  # I have simply set G = 1. In can in principle be changed.
  ###########################################################################
  if (length(FFQ_food_col) == 0) {
    with_covariates_ind <- 0
    nointake_ffq        <- numeric(0)
  }
  ###########################################################################
  # If you have an FFQ for the food, transform it. Also, generate a vector
  # that is the indicator that the FFQ for the food = 0
  ###########################################################################
  if (length(FFQ_food_col) > 0) {
    nointake_ffq            <- as.numeric(zz[ ,FFQ_food_col] == 0)
    nointake_ffq            <- matrix(nointake_ffq, nrow = n_subject, byrow = T)
    nointake_ffq            <- nointake_ffq[ ,1]
    zz[ ,FFQ_food_col]      <- boxcoxtrans1(zz[ ,FFQ_food_col], lambda_FFQ_food)
    ffq_food                <- zz[ ,FFQ_food_col]
    ffq_food                <- matrix(ffq_food, nrow = n_subject, byrow = T)
    ffq_food                <- ffq_food[ ,1]
  }
  
  ###########################################################################
  # If you have an FFQ for the food, transform it. The transformations were
  # saved when running All_consumers.m and do not need to be recalculated
  ###########################################################################
  if (length(FFQ_energy_col) > 0) {
    zz[ ,FFQ_energy_col]      <- boxcoxtrans1(zz[ ,FFQ_energy_col], lambda_FFQ_energy)
  }
  
  ###########################################################################
  # Now create the data to have one row per person
  ###########################################################################
  ID          <- zz[ ,ID_col]
  n_subject   <- length(unique(ID))
  ID_new      <- rep(1:n_subject,each = 4)
  for (jj      in 1:n_subject) {
    uu      <- (ID_new == jj)
    vv      <- zz[uu, ]
    if (jj  == 1) {
      xx  <- vv[1, X_cols]
      if (mmi == 2) {
        rec_food   <- c(vv[1,Food_col], vv[2,Food_col])
        rec_energy <- c(vv[1,Energy_col], vv[2,Energy_col])
      }
      if (mmi == 3) {
        rec_food   <- c(vv[1,Food_col], vv[2,Food_col], vv[3,Food_col])
        rec_energy <- c(vv[1,Energy_col], vv[2,Energy_col], vv[3,Energy_col])
      }
      if (mmi == 4) {
        rec_food   <- c(vv[1,Food_col], vv[2,Food_col], vv[3,Food_col], vv[4,Food_col])
        rec_energy <- c(vv[1,Energy_col], vv[2,Energy_col], vv[3,Energy_col], vv[4,Energy_col])
      }
    }
    if (jj  > 1) {
      xx  <- rbind(xx, vv[1,X_cols])
      if (mmi == 2) {
        rec_food   <- rbind(rec_food, c(vv[1,Food_col], vv[2,Food_col]))
        
        rec_energy <- rbind(rec_energy, c(vv[1,Energy_col], vv[2,Energy_col]))
        
      }
      if (mmi == 3) {
        rec_food   <- rbind(rec_food, c(vv[1,Food_col], vv[2,Food_col], vv[3,Food_col]))
        
        rec_energy <- rbind(rec_energy, c(vv[1,Energy_col], vv[2,Energy_col], vv[3,Energy_col]))
        
      }
      if (mmi == 4) {
        rec_food   <- rbind(rec_food, c(vv[1,Food_col], vv[2,Food_col], vv[3,Food_col], vv[4,Food_col]))
        
        rec_energy <- rbind(rec_energy, c(vv[1,Energy_col], vv[2,Energy_col], vv[3,Energy_col], vv[4,Energy_col]))
        
      }
    }
  }
  
  ###########################################################################
  # Standardize the covariate data. If there is an FFQ, it has already been
  # transformed
  ###########################################################################
  mm             <- ncol(xx)
  for (jj         in 1:mm) {
    mu         <- mean(xx[ ,jj])
    sig        <-  sd(xx[ ,jj])
    xx[ ,jj]   <- (xx[ ,jj]  - mu) / sig
  }
  ###########################################################################
  # Set up the design matrices.
  ###########################################################################
  n              <- nrow(xx)
  x_design       <- as.matrix(cbind(rep(1,n), xx))
  
  Xtildei <- array(1, c(n,ncol(x_design), 3))
  for (i in 1:3) {
    Xtildei[ , , i] = x_design
  }
  ###########################################################################
  # Initialization of beta_temp if its null
  ###########################################################################
  
  if (is.null(beta_temp)) {
    beta_temp = matrix(0, (mm + 1), 3)
  }
  beta_temp = round(beta_temp,4)
  
  consumer_percent_prior_mean <- ((with_covariates_ind == 0) * 0.8 + (with_covariates_ind != 0) * 0.5)
  # the prior mean of the percentage of consumers
  # 0.5 for with covariates in alpha;
  # 0.8 for without covariates in alpha
  consumer_percent_start <- ((with_covariates_ind == 0) * 0.8 + (with_covariates_ind != 0) * 0.5)
  # starting value for the percentage of consumers
  # 0.5 for with covariates in alpha;
  # 0.8 for without covariates in alpha
  alpha_prior_sd <- ((with_covariates_ind == 0) * 0.4 + (with_covariates_ind != 0) * 1)
  
  # prior standard deviation for alpha
  # 0.4 for with covariates in alpha;
  # 0.2 for without covariates in alpha
  myseed = 6309021
  set.seed(myseed)
  
  
  ###########################################################################
  # Do a simple processing of the data
  ###########################################################################
  pd <- process_data(rec_food, rec_energy, n, mmi, lambda_rec_food,
                     lambda_rec_energy)
  
  Wistar       <- pd$Wistar
  Wi2          <- round(pd$Wi2, 4)
  Wi3          <- round(pd$Wi3, 4)
  didconsume   <- pd$didconsume
  a0_food      <- pd$a0_food
  a0_energy    <- pd$a0_energy
  mumu         <- round(pd$mumu,4)
  sigsig       <- round(pd$sigsig,4)
  mu_e         <- round(pd$mu_e, 4)
  sig_e        <- round(pd$sig_e,4)
  
  ###########################################################################
  # Save half of the minimum positive food value, half of minimum energy
  # value, and Box-Cox transformation parameters.
  ###########################################################################
  # writeMat("NOutput_Females_4Recalls/a0.mat", a0_food = a0_food)
  # writeMat("NOutput_Females_4Recalls/a0_energy.mat", a0_energy = a0_energy)
  # writeMat("NOutput_Females_4Recalls/lambda_REC.mat", lambda_rec_food = lambda_rec_food)
  # writeMat("NOutput_Females_4Recalls/lambda_REC_Energy.mat", lambda_rec_energy = lambda_rec_energy)
  # ###########################################################################
  # Set up the design matrices in the evewr consumer model.
  ###########################################################################
  if (with_covariates_ind == 1) {
    GGalpha <- cbind(rep(1,n), ffq_food, nointake_ffq)
  }
  if (with_covariates_ind == 2) {
    GGalpha <- cbind(rep(1,n), ffq_food)
  }
  if (with_covariates_ind == 3) {
    GGalpha <- cbind(rep(1,n), nointake_ffq)
  }
  if (with_covariates_ind == 0) {
    GGalpha <- rep(1,n)
  }
  
  
  ###########################################################################
  ## Set the MCMC parameters
  ###########################################################################
  print('Set the MCMC parameters')
  ###########################################################################
  # Set the prior and starting value for alpha
  ###########################################################################
  prior_alpha_mean      <- qnorm(consumer_percent_prior_mean,0,1)*rep(1, ncol(GGalpha))
  
  prior_alpha_cov       <- (alpha_prior_sd ^ 2)*diag(rep(1, ncol(GGalpha)))
  alpha_start           <- qnorm(consumer_percent_start,0,1) * rep(1, ncol(GGalpha))
  
  alpha                 <- alpha_start
  
  ###########################################################################
  # Set the dimension of beta
  ###########################################################################
  mdesign  <- ncol(xx) + 1
  ###########################################################################
  # Set the prior and starting value for beta
  # Good starting values for the beta parameters and their covariance
  # matrices are available by other means. I ran the consumption program and
  # the amount program to get these values.
  ###########################################################################
  
  #print('check 1')
  if (beta_prior_mean_ind  == 9999) {
    prior_beta_mean  <- beta_temp
  } else if (beta_prior_mean_ind == 0) {
    prior_beta_mean  <- matrix(0,mdesign,3)
  }
  
  temp       <- (10 * diag(rep(1,mdesign)))
  prior_beta_cov <- array(1, c(nrow(temp),ncol(temp), 3))
  for (i in 1:3) {
    prior_beta_cov[ , , i] = temp
  }
  
  if (beta_start_ind == 9999) {
    beta_start <- beta_temp
  } else {
    beta_start <- beta_start_ind * rep(1, length(prior_beta_mean))
  }
  beta <- beta_start
  rm(beta_start)
  ###########################################################################
  # Set the prior and starting value for Sigmae
  ###########################################################################
  #print('check 2')
  r     <- 0
  theta <- 0
  s22   <- 1
  s33   <- 1
  R     <- matrix(c(1, 0, r*cos(theta),
                    0, 1,  r*sin(theta),
                    r*cos(theta),  r*sin(theta), 1), nrow = 3, byrow = T)
  
  A <- diag(c(1, sqrt(s22), sqrt(s33)))
  Sigmae       <- A %*% R %*% A
  iSigmae      <- solve(Sigmae)
  Sigmae_start <- Sigmae
  prior_Sigmae_doff <- 5
  ###########################################################################
  # Set the prior and starting values for Sigmau
  ###########################################################################
  #print('check 3')
  prior_Sigmau_mean  <- diag(3)
  sig    <- aarp_setprior_Sigmau_never(prior_Sigmau_mean)
  Sigmau_temp_regular <- sig$prior_Sigmau_mean
  prior_Sigmau_doff   <- sig$prior_Sigmau_doff
  
  if (Sigmau_prior_mean_ind == 1) {
    prior_Sigmau_mean <- Sigmau_temp_episodically
  }else if (Sigmau_prior_mean_ind == 2) {
    prior_Sigmau_mean <- Sigmau_temp_regular
  }else if (Sigmau_prior_mean_ind == 3) {
    prior_Sigmau_mean <- Sigmau_temp_regular / 2
  }else{
    stop(paste("Sigmau_prior_mean_ind not recognized"))
  }
  
  if (Sigmau_start_ind == 1) {
    Sigmau_start <- Sigmau_temp_episodically
  }else if (Sigmau_start_ind == 2) {
    Sigmau_start <- Sigmau_temp_regular
  }else{
    stop(paste("Sigmau_start_ind not recognized"))
  }
  
  Sigmau    <- round(Sigmau_start, 4)
  iSigmau    <- round(pracma::inv(Sigmau), 4)
  
  ###########################################################################
  ## Initialize a few things
  ###########################################################################
  ###########################################################################
  # Set starting values for the Utildei
  ###########################################################################
  #set.seed(71094)
  #print('check 4')
  Utildei <- pracma::randn(n,ncol(Sigmau)) %*% pracma::sqrtm(Sigmau)$B
  ###########################################################################
  # Get starting values for the W_{ijk}
  ###########################################################################
  
  WtildeiS = array(1, c(n, 3, 4))
  WtildeiS[ ,2, ]  <- Wi2
  WtildeiS[ ,3, ]  <- Wi3
  
  WtildeiS[ ,1, ]  <- abs(pracma::repmat(Xtildei[ , ,1] %*% beta[ ,1] + Utildei[ ,1], 1,mmi)
                          + pracma::randn(n,mmi))
  WtildeiS[ ,1, ]  <- (WtildeiS[ ,1, ] * Wistar) - (WtildeiS[ ,1, ] * (1 - Wistar))
  Wtildei          <- WtildeiS
  numgen           <- 20
  Wtildeinew       <- gen_Wtildei_1foodplusenergy_never(WtildeiS,beta,Xtildei,
                                                        Utildei,n,iSigmae,Wistar,mmi, numgen)
  Wtildei          <- round(Wtildeinew, 4)
  Wtildei_start <- Wtildei
  #print('check 5')
  
  ###########################################################################
  ## MCMC
  ###########################################################################
  
  start.time <- Sys.time()
  
  MCMC_analysis = perform_MCMC(nMCMC,nthin, n, mmi, Xtildei, Utildei, beta, alpha, GGalpha, didconsume, 
                               prior_alpha_cov, prior_alpha_mean, Wtildei, iSigmae, iSigmau, Wistar,r, theta,
                               s22, s33, Sigmau, Sigmae, prior_Sigmau_mean, prior_Sigmau_doff, prior_Sigmae_doff,
                               prior_beta_mean, prior_beta_cov, update_beta1_var_ind, lambda_rec_food, 
                               lambda_rec_energy, mumu, sigsig, mu_e, sig_e, mdesign, rw_ind, beta1_accept_count, 
                               a0_food, a0_energy, ndist, ndim)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
  # 
  # library(microbenchmark)
  # microbenchmark(perform_MCMC(nMCMC,nthin, n, mmi, Xtildei, Utildei, beta, alpha, GGalpha, didconsume, 
  #                             prior_alpha_cov, prior_alpha_mean, Wtildei, iSigmae, iSigmau, Wistar,r, theta,
  #                             s22, s33, Sigmau, Sigmae, prior_Sigmau_mean, prior_Sigmau_doff, prior_Sigmae_doff,
  #                             prior_beta_mean, prior_beta_cov, update_beta1_var_ind, lambda_rec_food, 
  #                             lambda_rec_energy, mumu, sigsig, mu_e, sig_e, mdesign, rw_ind, beta1_accept_count, 
  #                             a0_food, a0_energy, ndist, ndim), times = 10)
  # # median of 53 mins
  
  
  alpha_trace = MCMC_analysis$alpha_trace
  beta_trace = MCMC_analysis$beta_trace
  never_trace = MCMC_analysis$never_trace
  r_trace = MCMC_analysis$r_trace
  theta_trace = MCMC_analysis$theta_trace
  s22_trace = MCMC_analysis$s22_trace
  s33_trace = MCMC_analysis$s33_trace
  Sigmae_trace = MCMC_analysis$Sigmae_trace
  Sigmau_trace = MCMC_analysis$Sigmau_trace
  usual_intake_energy_trace = MCMC_analysis$usual_intake_energy_trace
  usual_intake_food_trace = MCMC_analysis$usual_intake_food_trace
  ###########################################################################
  ## Thinning, burn-in, compute posterior mean, standard deviation,
  # credible interval, Compute distribution of usual intake, save results
  ###########################################################################
  ###########################################################################
  # Thinning and save the traces (after thinning, but before burn-in)
  ###########################################################################
  
  nn            <- nrow(r_trace)
  mm            <- floor((nn + eps) / nthin)
  thin_index    <- seq(nthin,(mm*nthin), by = nthin)
  alpha_thin_trace   <-  alpha_trace[thin_index, ]
  never_thin_trace   <-  never_trace[thin_index,1]
  r_thin_trace   <-      r_trace[thin_index,1]
  theta_thin_trace   <-  theta_trace[thin_index,1]
  s22_thin_trace   <-    s22_trace[thin_index,1]
  s33_thin_trace   <-    s33_trace[thin_index,1]
  Sigmae_thin_trace <- Sigmae_trace[ , ,thin_index]
  Sigmau_thin_trace <- Sigmau_trace[ , ,thin_index]
  beta_thin_trace <-   beta_trace[ , ,thin_index]
  
  
  ###########################################################################
  # Get rid of the burn-in
  ###########################################################################
  mmburn             <- floor((nburn + eps) / nthin)
  alpha_thin_trace   <-  alpha_thin_trace[(mmburn + 1):mm, ]
  never_thin_trace   <-  never_thin_trace[(mmburn + 1):mm]
  r_thin_trace   <-      r_thin_trace[(mmburn + 1):mm]
  theta_thin_trace   <-  theta_thin_trace[(mmburn + 1):mm]
  s22_thin_trace   <-    s22_thin_trace[(mmburn + 1):mm]
  s33_thin_trace   <-    s33_thin_trace[(mmburn + 1):mm]
  Sigmae_thin_trace   <- Sigmae_thin_trace[ , ,(mmburn + 1):mm]
  Sigmau_thin_trace   <- Sigmau_thin_trace[ , ,(mmburn + 1):mm]
  beta_thin_trace   <-   beta_thin_trace[ , ,(mmburn + 1):mm]
  ###########################################################################
  # Compute and save the posterior means, standard deviation and
  # 95# credible interval
  ###########################################################################
  alpha_postmean      <-  apply(alpha_thin_trace, 2, mean) #colMeans(alpha_thin_trace)
  never_postmean      <-  mean(never_thin_trace)
  beta_postmean       <-  apply(beta_thin_trace, c(1,2), mean)
  Sigmau_postmean     <-  apply(Sigmau_thin_trace, c(1,2), mean)
  Sigmae_postmean     <-  apply(Sigmae_thin_trace, c(1,2), mean)
  
  alpha_postsd        <- apply(alpha_thin_trace, 2, sd)
  never_postsd        <- sd(never_thin_trace)
  beta_postsd         <- apply(beta_thin_trace, c(1,2), sd)
  Sigmau_postsd       <- apply(Sigmau_thin_trace, c(1,2), sd)
  Sigmae_postsd       <- apply(Sigmae_thin_trace, c(1,2), sd)
  
  alpha_ci            <- apply(alpha_thin_trace, 2, function(x) quantile(x, c(0.025, 0.975)))
  never_ci            <- quantile(never_thin_trace,c(0.025, 0.975))
  beta_ci             <- apply(beta_thin_trace, c(1,2), function(x) quantile(x, c(0.025, 0.975)))
  Sigmau_ci           <- apply(Sigmau_thin_trace, c(1,2), function(x) quantile(x, c(0.025, 0.975)))
  Sigmae_ci           <- apply(Sigmae_thin_trace, c(1,2), function(x) quantile(x, c(0.025, 0.975)))
  
  # Computer the correlation matrix for Sigmau
  Corru_postmean <- matrix(NaN, nrow(Sigmau_postmean), ncol(Sigmau_postmean))
  for (iicorr in 1:nrow(Sigmau_postmean)) {
    for (jjcorr in 1:ncol(Sigmau_postmean)) {
      Corru_postmean[iicorr,jjcorr] <- Sigmau_postmean[iicorr,jjcorr] /
        sqrt(Sigmau_postmean[iicorr,iicorr]*Sigmau_postmean[jjcorr,jjcorr])
    }
  }
  # Computer the correlation matrix for Sigmae
  Corre_postmean <- matrix(NaN, nrow(Sigmae_postmean), ncol(Sigmae_postmean))
  for (iicorr in 1:nrow(Sigmae_postmean)) {
    for (jjcorr in 1:ncol(Sigmae_postmean)) {
      Corre_postmean[iicorr,jjcorr] <- Sigmae_postmean[iicorr,jjcorr] /
        sqrt(Sigmae_postmean[iicorr,iicorr]*Sigmae_postmean[jjcorr,jjcorr])
    }
  }
  
  ###########################################################################
  # Compute distribution of usual intake of food, energy and
  # food/(energy/1000).
  # Compute the cdf of usual intake for consumers on a fine grid.
  # We have an estimate of the cdf for consumers,
  # which can be inverted to get the percentiles.
  ###########################################################################
  usual_intake_ratio_trace <- 1000 * usual_intake_food_trace /
    usual_intake_energy_trace

  #food_p_mat = [0:0.0001:0.2, 0.2005:0.0005:1]';
  food_p_mat <- seq(0,1,length.out = 501)
  aa <- usual_intake_food_trace[!is.nan(usual_intake_food_trace)]
  food_distribution <- make_percentiles_without_weight(aa, food_p_mat)
  mu_ui_food <- mean(aa)
  sig_ui_food <- sd(aa)
  ##plot(seq(0.01, 0.01, 0.99), food_distribution)

  # #energy_p_mat = [0:0.0001:1]';
  energy_p_mat <- seq(0,1,length.out = 501)
  aa <- usual_intake_energy_trace[!is.nan(usual_intake_energy_trace)]
  #aa <- (aa-min(aa))/(max(aa)-min(aa))
  energy_distribution <- make_percentiles_without_weight(aa, energy_p_mat)
  mu_ui_energy = mean(aa)
  sig_ui_energy = sd(aa)
  # #plot(0.01:0.01:0.99, energy_distribution)

  #ratio_p_mat = [0:0.0001:0.2, 0.2005:0.0005:1]';
  ratio_p_mat = seq(0,1,length.out = 501)
  aa = usual_intake_ratio_trace[!is.nan(usual_intake_ratio_trace)]
  ratio_distribution = make_percentiles_without_weight(aa, ratio_p_mat)
  mu_ui_ratio = mean(aa)
  sig_ui_ratio = sd(aa)
  # # #plot(0.01:0.01:0.9mu_ui_energy <- mean(aa)

  adj_factor = 1  # This is in here because the Norfolk energy data are on a
  # different scale from the EATS energy data.
  ui_percentile_ind = c(5, 10, 25, 50, 75, 90, 95)
  food_distribution_percentile = t(food_distribution[ui_percentile_ind]);
  energy_distribution_percentile = t(energy_distribution[ui_percentile_ind])/adj_factor;
  ratio_distribution_percentile = t(ratio_distribution[ui_percentile_ind])*adj_factor;

  # ###########################################################################
  # # What percentage of MCMC steps in which beta1 moves
  # ###########################################################################
  beta1_accept_rate <- beta1_accept_count/nMCMC
  # disp(paste('beta1 accept rate is ', num2str(beta1_accept_rate*100), '#'))
  
  
  print(' ')
  print(' ')
  print(' ')
  print(' ')
  print(' ')
  print(' ')
  print(' ')
  print(' ')
  
  
  beta_11 <-  beta_thin_trace[ 1,1, ]
  beta_12 <-  beta_thin_trace[ 2,1, ]
  beta_13 <-  beta_thin_trace[ 3,1, ]
  beta_14 <-  beta_thin_trace[ 4,1, ]
  beta_15 <-  beta_thin_trace[ 5,1, ]
  beta_21 <-  beta_thin_trace[ 1,2, ]
  beta_22 <-  beta_thin_trace[ 2,2, ]
  beta_23 <-  beta_thin_trace[ 3,2, ]
  beta_24 <-  beta_thin_trace[ 4,2, ]
  beta_25 <-  beta_thin_trace[ 5,2, ]
  beta_31 <-  beta_thin_trace[ 1,3, ]
  beta_32 <-  beta_thin_trace[ 2,3, ]
  beta_33 <-  beta_thin_trace[ 3,3, ]
  beta_34 <-  beta_thin_trace[ 4,3, ]
  beta_35 <-  beta_thin_trace[ 5,3, ]
  
  par(mfrow = c(1,1))
  matplot(t(beta_thin_trace[ ,1, ]), type = "l", ylab = expression(beta[1]))
  nn <- ncol(t(beta_thin_trace[ ,1, ]))
  legend("topleft", legend = c(expression(beta[11]), expression(beta[12]), expression(beta[13]), expression(beta[14]), expression(beta[15])),col = seq_len(nn),cex = 0.8,fill = seq_len(nn))
  title(expression(paste("Traceplot of ", beta[1])))

  matplot(t(beta_thin_trace[ ,2, ]), type = "l", ylab = expression(beta[2]))
  nn <- ncol(t(beta_thin_trace[ ,2, ]))
  legend("topleft", legend = c(expression(beta[21]), expression(beta[22]), expression(beta[23]), expression(beta[24]), expression(beta[25])),col = seq_len(nn),cex = 0.8,fill = seq_len(nn))
  title(expression(paste("Traceplot of ",beta[2])))

  matplot(t(beta_thin_trace[ ,3, ]), type = "l", ylab = expression(beta[3]))
  nn <- ncol(t(beta_thin_trace[ ,3, ]))
  legend("topleft", legend = c(expression(beta[31]), expression(beta[32]), expression(beta[33]), expression(beta[34]), expression(beta[35])),col = seq_len(nn),cex = 0.8,fill = seq_len(nn))
  title(expression(paste("Traceplot of ",beta[3])))

  par(mfrow = c(2,1))
  plot(alpha_thin_trace[ ,1], type = "l", ylab = expression(alpha[1]))
  title(expression(paste("Traceplot of ",alpha[1])))
  plot(alpha_thin_trace[ ,2], type = "l", ylab = expression(alpha[2]))
  title(expression(paste("Traceplot of ",alpha[2])))

  par(mfrow = c(1,1))
  plot(never_thin_trace, type = "l", ylab = "Probability")
  title("Traceplot of probability of being never consumers")

  par(mfrow = c(2,2))
  plot(r_thin_trace, type = "l", ylab = "r")
  title(expression(paste("Traceplot of ",r)))
  plot(theta_thin_trace, type = "l", ylab = expression(theta))
  title(expression(paste("Traceplot of ",theta)))
  plot(s22_thin_trace, type = "l", ylab = expression(s[22]))
  title(expression(paste("Traceplot of ",s[22])))
  plot(s33_thin_trace, type = "l", ylab = expression(s[33]))
  title(expression(paste("Traceplot of ",s[33])))


  par(mfrow = c(3,3))
  plot((Sigmae_thin_trace[1,1,]), type = "l", ylab = expression(sigma[e[11]]),xlab = "Iteration")
  title(expression(paste("Traceplot of ",sigma[e[11]])))
  plot((Sigmae_thin_trace[1,2,]), type = "l", ylab = expression(sigma[e[12]]),xlab = "Iteration")
  title(expression(paste("Traceplot of ",sigma[e[12]])))
  plot((Sigmae_thin_trace[1,3,]), type = "l", ylab = expression(sigma[e[13]]),xlab = "Iteration")
  title(expression(paste("Traceplot of ",sigma[e[13]])))
  plot((Sigmae_thin_trace[2,1,]), type = "l", ylab = expression(sigma[e[21]]),xlab = "Iteration")
  title(expression(paste("Traceplot of ",sigma[e[21]])))
  plot((Sigmae_thin_trace[2,2,]), type = "l", ylab = expression(sigma[e[22]]),xlab = "Iteration")
  title(expression(paste("Traceplot of ",sigma[e[22]])))
  plot((Sigmae_thin_trace[2,3,]), type = "l", ylab = expression(sigma[e[23]]),xlab = "Iteration")
  title(expression(paste("Traceplot of ",sigma[e[23]])))
  plot((Sigmae_thin_trace[3,1,]), type = "l", ylab = expression(sigma[e[31]]),xlab = "Iteration")
  title(expression(paste("Traceplot of ",sigma[e[31]])))
  plot((Sigmae_thin_trace[3,2,]), type = "l", ylab = expression(sigma[e[32]]),xlab = "Iteration")
  title(expression(paste("Traceplot of ",sigma[e[32]])))
  plot((Sigmae_thin_trace[3,3,]), type = "l", ylab = expression(sigma[u[33]]),xlab = "Iteration")
  title(expression(paste("Traceplot of ",sigma[e[33]])))

  par(mfrow = c(3,3))
  plot((Sigmau_thin_trace[1,1,]), type = "l", ylab = expression(sigma[u[11]]),xlab = "Iteration")
  title(expression(paste("Traceplot of ",sigma[u[11]])))
  plot((Sigmau_thin_trace[1,2,]), type = "l", ylab = expression(sigma[u[12]]),xlab = "Iteration")
  title(expression(paste("Traceplot of ",sigma[u[12]])))
  plot((Sigmau_thin_trace[1,3,]), type = "l", ylab = expression(sigma[u[13]]),xlab = "Iteration")
  title(expression(paste("Traceplot of ",sigma[u[13]])))
  plot((Sigmau_thin_trace[2,1,]), type = "l", ylab = expression(sigma[u[21]]),xlab = "Iteration")
  title(expression(paste("Traceplot of ",sigma[u[21]])))
  plot((Sigmau_thin_trace[2,2,]), type = "l", ylab = expression(sigma[u[22]]),xlab = "Iteration")
  title(expression(paste("Traceplot of ",sigma[u[22]])))
  plot((Sigmau_thin_trace[2,3,]), type = "l", ylab = expression(sigma[u[23]]),xlab = "Iteration")
  title(expression(paste("Traceplot of ",sigma[u[23]])))
  plot((Sigmau_thin_trace[3,1,]), type = "l", ylab = expression(sigma[u[31]]),xlab = "Iteration")
  title(expression(paste("Traceplot of ",sigma[u[31]])))
  plot((Sigmau_thin_trace[3,2,]), type = "l", ylab = expression(sigma[u[32]]),xlab = "Iteration")
  title(expression(paste("Traceplot of ",sigma[u[32]])))
  plot((Sigmau_thin_trace[3,3,]), type = "l", ylab = expression(sigma[u[33]]),xlab = "Iteration")
  title(expression(paste("Traceplot of ",sigma[u[33]])))

  
  
  return(list(alpha_postmean = alpha_postmean,alpha_postsd = alpha_postsd, alpha_ci = alpha_ci, 
              never_postmean = never_postmean, never_postsd = never_postsd, never_ci = never_ci,
              beta_postmean = beta_postmean, beta_postsd = beta_postsd, beta_ci = beta_ci, 
              Sigmau_postmean = Sigmau_postmean, Sigmau_postsd = Sigmau_postsd, Sigmau_ci = Sigmau_ci, 
              Sigmae_postmean = Sigmae_postmean, Sigmae_postsd = Sigmae_postsd, Sigmae_ci = Sigmae_ci, 
              mu_ui_food = mu_ui_food, sig_ui_food = sig_ui_food, mu_ui_energy = mu_ui_energy,
              sig_ui_energy = sig_ui_energy, mu_ui_ratio = mu_ui_ratio, sig_ui_ratio = sig_ui_ratio,
              food_distribution_percentile = food_distribution_percentile, energy_distribution_percentile = energy_distribution_percentile,
              ratio_distribution_percentile = ratio_distribution_percentile))
}

