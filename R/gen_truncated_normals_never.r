gen_truncated_normals_never <- function(trunc_value,startxi,numgen){
#
# This generates standard normals truncated from the left
# at trunc_value, with starting values startxi,
# using a rejection sampler devised by C. P. Robert
# Statistics and Computing, 1995, volumn 5, pp 121-125.
#
# The rejection sampler is only used for truncation values
# > 0, because it is remarkably bad for truncation values
# < 0, i.e., it tends to stay where it starts in this case.
#
# INPUT:
#      trunc_value: the random variable is truncated
#                   at the left from trunc_value, a vector
#      startxi:     starting values
#      numgen:      number of times you try (recommended = 50)
#
# OUTPUT:
#      truncated normals, the same dimension as trunc_value
#
#
# library(iemisc)
eps    <- 2.2204e-16
n      <- nrow(trunc_value)
alpha  <- (trunc_value + sqrt(4 + (trunc_value ^ 2))) / 2
thesign <- as.numeric(trunc_value >= 0)# Whether the truncation point
                              # is positive
genww  <- trunc_value * as.numeric(trunc_value > 0)
temp2  <- pracma::randn(n,1)
for (jj in 1:numgen){
    xicand <- trunc_value - ( (1 / alpha) * log(pracma::rand(n,1)))
    mmmm   <- (pracma::randn(n,1) < exp(-.5 * ( (xicand - alpha) ^ 2)))
    temp1  <- (xicand  * as.numeric(mmmm == 1)) + (genww * as.numeric(mmmm == 0))
    ssss   <- runif(n,-1,1)
    temp2  <- (temp2 * as.numeric(ssss < trunc_value)) +
               (ssss * as.numeric(ssss >= trunc_value))
    genww  <- (temp2 * as.numeric(thesign == 0)) + (temp1 *
                                                    as.numeric(thesign == 1))
}
genww  <- (genww * (genww > trunc_value)) + ((trunc_value + eps)
                             *(genww <= trunc_value))
return(genww)
}



