backtransform<- function(lambda,Xtildei_1,beta_1,sigmae_1,mumu,sigsig,Utildei_1,n){
# Compute the 9 point backtransformation for any component
#
# INPUT
# lambda      = the tranformation parameter
# Xtildei     = the design matrix
# beta        = the posterior mean
# sigmae      = standard deviation of epsilon
# mumu        = the transformed mean
# sigsig      = the transformed standard deviation
# Utildei     = the realized U
# n           = the sample size

# set the abscissas and weights
x <- c(-2.1,-1.3,-0.8,-0.5, 0.00, 0.5, 0.8, 1.3, 2.1)
w <- c(0.063345, 0.080255, 0.070458, 0.159698, 0.252489,0.159698,
       0.070458, 0.080255, 0.063345)

    temp <- (pracma::repmat(Xtildei_1 %*% beta_1 + as.matrix(Utildei_1), 1,pracma::size(x,2))
        +  pracma::repmat(x, n, 1) * sqrt(2) * sigmae_1)# get the terms in the innermost parentheses
    temp <-  (mumu + sigsig * temp / sqrt(2))# de-standardize
    temp <- MASS::ginv(temp,lambda)

    temp <- colSums(t(pracma::repmat(w,n,1))*temp)#get n*1 vector,  eq A.5

    return(temp)

 }
