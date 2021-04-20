#function s33new = updated_parameter_s33(r,theta,s22,s33,Xtildei,Utildei,Wtildei,n,beta);
updated_parameter_s33_never <- function(r,theta,s22,s33,qq,mmi,n){

# Do the Metropolis Step for s33.
# r = current value
# theta = current value
# s22 = current value
# s33 = current value

s33curr       <- s33
s33cand       <- s33 + (0.4 * (runif(1) - 0.5))
GofSigmaecurr = formGofSigmae_never_c(r,theta,s22,s33curr,qq,mmi)
GofSigmaecand = formGofSigmae_never_c(r,theta,s22,s33cand,qq,mmi)

gg            <- GofSigmaecand - GofSigmaecurr
gg            <- gg - ((mmi*n/2) * log(s33cand)) + ((mmi*n/2) * log(s33curr))
gg            <- min(1, exp(gg)*as.numeric((s33cand>=0)&(s33cand<=3)))
ss            <- pracma::rand(1,1)
s33new       <- (s33cand * (ss < gg)) +  (s33curr * (ss > gg))

return(s33new)

}

