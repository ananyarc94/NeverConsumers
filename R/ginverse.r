
ginverse <- function(z,lambda){
#
# Computes the inverse of the Box-Cox transformation, and makes sure that
# the argument does not get negative
#
lambda = as.numeric(lambda)  
if (lambda == 0){
    x <- exp(z)
}
if (lambda > 0){
    w <- 1 + (lambda * z)
    w[w<0] <- 0
    x <- w ^ (1 / lambda)
}

return(x)

}
