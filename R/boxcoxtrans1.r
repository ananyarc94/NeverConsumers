
boxcoxtrans1 <- function(x, lam){
###########################################################################
# This function computes Box-Cox transformation described in
# Transformations in the AARP.
# input:    x:          vector for Box-Cox transformatioin
#           lam:        prespecified lambda, if you know it already, otherwise it is estimated by the function
# output:   bcx:        x after Box-Cox transformation
#           bclambda:   lambda that maximize R^2
###########################################################################
lam = as.numeric(lam)

if (length(x) == 1){
    warning('x must be a vector.')
}
if (any(x < 0)){
    warning('x must be nonnegative.')
}

    nargin <- nargs()
###########################################################################
# if only one input, then we need to estimate lambda
###########################################################################
if (nargin == 1){
    ind = which(x > 0)
    x <- x[ind] # take the non-zero data
    #n = size(x,1);
    n <- 99 # corrected by Nelis
    lambda <- seq(0, 1, 0.01)
    nlambda <- length(lambda)
    p <- seq(0.01, 0.99, 0.01)
    xq <- quantile(x, p)
    # i = linspace(1,n,101);
    # i = i(2:(end-1));
    i <- 1:99
    b <- qnorm((i-3/8)/(n+1/4))
    rsq <- numeric(nlambda)

    for (j in 1:nlambda){
        if (lambda[j] == 0){
            y <- log(xq)
        } else {
            y <- (xq^(lambda[j])-1)/lambda[j]
        }

        rsq[j] <- summary(lm(y ~ b))$adj.r.squared
    }


    bclambda <- lambda[which.max(rsq)]

    #plot(lambda,rsq)

###########################################################################
# if lambda is prespecifed, then only do the transformation.
###########################################################################
} else if (nargin == 2){
    if (length(lam) != 1){
        stop(paste('lambda must be a scalar.'))
    }
    bclambda <- lam
    #else error('Only up to two inputs are supported.')
    if (any(any(x<=0)) && (bclambda == 0)){
        a0 <- min(min(x[x>0]))/2
        x[x<=0] <- a0
        print(paste('There are non-positive x and the lambda is ', (bclambda)))
        print(paste('To solve the problem, the non-positive values are replaced by '))
        print(paste('half of minimum of the positive values, which is ',(a0),'.'))
     }
}

if (bclambda == 0){
    bcx <- log(x)
} else {
    bcx <- (x^(bclambda)-1)/bclambda
}

    return(bcx)
}
