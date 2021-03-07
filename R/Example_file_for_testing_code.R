###########################################################################
# This is the MCMC analysis for a single food plus energy in the three-part
# model. 
# Last revised by Ananya Roy Chowdhury <ananya@stat.tamu.edu> on 02/28/2021.
# MATLAB Version: 8.1.0.604 (R2013a)
###########################################################################

# These data have been cleaned. The way this is done is that any person who
# has an "outlier" in the recalls for FFQs is deleted from the analysis
# (this is because the program requires the same number of recalls per
# person). An outlier is defined as follows. First transform the data using
# the program boxcoxtrans_20130925.m. Then, for the food, an outlier is any
# measurement that is more than two IQR about the 75th percentile. For
# energy, the same thing, but also any person who is more than 2 IQR below
# the 25th percentile.
# 
###########################################################################
## Initial Setup
###########################################################################
rm(list = ls())

setwd("C:\\Users\\anany\\Desktop\\NeverConsumers_NEW\\Supplemental_Material_Programs")
library(icesTAF)
# 
# source('boxcoxtrans1.r')
# source('process_data_r.r')
# source('aarp_setprior_Sigmau_never.r')
# source('gen_Wtildei_1foodplusenergy_never.r')
# source('gen_truncated_normals_never.r')
# source('update_Ni_with_covariates.r')
# source('gen_Wtildei_1foodplusenergy_never.r')
# source('updated_parameter_r_never.r')
# source('updated_parameter_theta_never.r')
# source('updated_parameter_s22_never.r')
# source('updated_parameter_s33_never.r')
# source('update_iSigmau.r')
# source('update_Utildei.r')
# source('update_beta1_with_prior_mean_random_walk.r')
# source('update_beta1_with_prior_mean.r')
# source('update_beta2_with_prior_mean.r')
# source('update_beta3_with_prior_mean.r')
# source('perform_MCMC.r')
# source('gen_truncated_normals_never.r')
# source('formGofSigmae_never.r')
# source('backtransform.r')
# source('make_percentiles_without_weight.r')
# source('neverconsumers_1.r')

mmi             = 4          # Number of recalls, integer 
nMCMC           = 50000     # Number of MCMC iterations. The more the better
nburn           = 10000      # Size of the burn-in
nthin           = 50         # thining
ndist           = 200        # average of the last ndist MCMC steps (after thinning)
# to get the cdf of usual intake 

data_set = read.csv('Example_data\\Simulated_Females_4Recalls.csv', header = T) # Name of the data set you will load 
output_name   = 'NOutput_Females_4Recalls' # Name of the directory where the output goes

# from the All consumer model resides
distribution_name = 'Output_UsualIntake_4recalls_Females.txt'
# Name of the Latex file that gives various
# percentiles
summary_name = 'Output_MCMC_4recalls_Females.txt'
# Name of the Latex file that gives the MCMC output



X_cols          = 2:5      # Columns where the covariates are
Food_col        = 6        # Column where the episodic variable is
Energy_col      = 7        # Column where energy (continuous) is
ID_col          = 1        # Column for ID
FFQ_food_col    = 4        # Column where the FFQ for the food is. Set = [] if no FFQ
# Please be aware that if there is a FFQ, you
# really should add this in, because FFQ have
# a tendency to generate massive, high
# leverage outliers, as in the EATS data
FFQ_energy_col  = 5        # Column where the FFQ for energy is. Set = [] if no FFQ
# Please be aware that if there is a FFQ, you
# really should add this in, because FFQ have
# a tendency to generate massive, high
# leverage outliers, as in the EATS data
with_covariates_ind = 3      # What to include as covariates in the
# ever consumer model 
# 0: a column of ones. 
# 1: a column of ones, the FFQ, and 
#    the indicator that the FFQ=0. 
# 2: a column of ones and the FFQ.
# 3: a column of ones and the indicator 
#    that the FFQ=0. 
n_gen = 30                    # Number of realizations of usual intake to
# generate. Must be a positive integer, 
# or 0 if no realizations to be generated. 
if(ndist * nthin > nMCMC){
  print('Please decrease ndist or increase nMCMC')
}


mkdir(output_name)



###########################################################################
# Read in the transformation parameters from the All_consumer.m
###########################################################################
lambda_rec_food   = R.matlab::readMat("All_consumers/Output_4Recalls/lambda_rec_food.mat")
lambda_rec_food   = lambda_rec_food$lambda.rec.food
lambda_rec_food   = as.numeric(lambda_rec_food)

lambda_rec_energy = R.matlab::readMat("All_consumers/Output_4Recalls/lambda_rec_energy.mat") 
lambda_rec_energy   = lambda_rec_energy$lambda.rec.energy
lambda_rec_energy   = as.numeric(lambda_rec_energy)


if (length(FFQ_food_col) > 0){
  lambda_FFQ_food   = R.matlab::readMat("All_consumers/Output_4Recalls/lambda_FFQ_food.mat")
  lambda_FFQ_food   = lambda_FFQ_food$lambda.FFQ.food
  lambda_FFQ_food   = as.numeric(lambda_FFQ_food)
  
}
if (length(FFQ_energy_col) > 0){
  lambda_FFQ_energy = R.matlab::readMat("All_consumers/Output_4Recalls/lambda_FFQ_energy.mat") 
  lambda_FFQ_energy   = lambda_FFQ_energy$lambda.FFQ.energy
  lambda_FFQ_energy   = as.numeric(lambda_FFQ_energy)
  
}

#For initial values of beta and Sigmau, they are matrices which we read from the attached file
beta_temp <- R.matlab::readMat("All_consumers/Output_4Recalls/beta_postmean.mat")
beta_temp <- beta_temp$beta.postmean

Sigmau_temp_episodically <-R.matlab::readMat("All_consumers/Output_4Recalls/Sigmau_postmean.mat")


Sigmau_temp_episodically <- Sigmau_temp_episodically$Sigmau.postmean

beta_start_ind = 9999      # the starting value of beta
# 9999: episodically; 
# all other number: a 5*3 matrix with the 
# number as each cell
beta_prior_mean_ind  = 0   # the prior mean for beta
# 9999: episodically (estimated by Saijuan's code)
                                # 0: regular (a 5*3 matrix with all elements = 0) 
rw_ind = 1                 # do you want to use the random walk proposal
                            # for beta_1
                                # 1: yes, use randon walk proposal
                                # Normal(\beta_{1,\curr}, \C_2 /M)
                                # 0: no, use Normal(\C_2 \C_1, \C_2 /M) 
update_beta1_var_ind = 0.5 # the variance for updating beta1, i.e. C1 
                            # in section A.9, this is the M in 
                            # Normal(\beta_{1,\curr}, \C_2 /M)
                            # and Normal(\C_2 \C_1, \C_2 /M) 
Sigmau_start_ind = 1       # the starting value of Sigmau
                                # 1: episodically (estimated by Saijuan's code)
# 2: regular (a 3*3 matrix with diagonal 
              # elements = 1 and off-diagonal elements = 0.5) 
Sigmau_prior_mean_ind = 2  # the prior mean for Sigmau 
# 1: episodically (estimated by Saijuan's code)
                                # 2: regular (a 3*3 matrix with diagonal 
                                # elements = 1 and off-diagonal elements = 0.5) 
                                # 3: half of regular
ndim     = 3 # the number of dimensions, here = 3 for indicator, amount 
            # and energy
beta1_accept_count = 0 # count how many times beta1 moves

eps = 10^-4

nc_object = neverConsumers(data_set, mmi, X_cols, Food_col, Energy_col, ID_col, FFQ_food_col,
                          FFQ_energy_col, with_covariates_ind, n_gen, lambda_rec_food, 
                          lambda_rec_energy, lambda_FFQ_food, lambda_FFQ_energy,beta_temp,
                          Sigmau_temp_episodically, nburn, nMCMC, nthin, ndist, 
                          beta_start_ind, beta_prior_mean_ind, rw_ind, update_beta1_var_ind, 
                          Sigmau_start_ind, Sigmau_prior_mean_ind, ndim,beta1_accept_count, 
                          myseed,eps)

