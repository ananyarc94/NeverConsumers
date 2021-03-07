make_percentiles_without_weight <- function(xx, pgrid){
###########################################################################
# This program takes an N by B matrix xx and computes it first through 99th
# percentile, weighted by weights. The program is slow, but very steady
#
# INPUT
# xx       = an n x B matrix of data. Each column has the "bootstrap" data
#
# OUTPUT
# percentiles = the estimated percentiles using sampling weights
###########################################################################

# Find the minimum and maximum values of xx
minmin <- min(xx)
maxmax <- max(xx)
# Create a grid upon which to compute the empirical cdf
# ngrid  = 501;
ngrid <- length(pgrid)
ygrid  <- numeric(ngrid)
# xgrid  = linspace(minmin,maxmax,ngrid)'; #'
xgrid <- as.vector(minmin + (maxmax - minmin) * pgrid)
# Now reshape the n x B matrix xx. You stack the second column below the
# first column, the thirs below the second, etc.
#aa     = reshape(xx,n.*B,1);
# Now make the weights compatible with the reshaped data. This means you
# simply stack the set on n weights on itself B times.
# ww     = kron(ones(B,1),weights);
# You now have weights and data. Compute the weighted empirical cdf on the
# grid
for(jj in 1:ngrid){
    ygrid[jj] <- sum(xx <= (xgrid[jj])) / length(xx)
}
# Now you have the weighted empirical cdf. The next step is to work out the
# percentiles. You are looking for the value of x grid that comes closest
# to giving any percentile
pp <- t(pracma::linspace(1,99,99))/ 100#'

vv <- 0.0001 * pracma::rand(ngrid,1)
vv <- pracma::sortrows(vv,1)

percentiles <- pracma::interp1(c(0,ygrid+vv),c(0,xgrid),c(0,pp),'cubic')

return(percentiles)

}
