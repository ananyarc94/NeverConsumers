#function s22new = updated_parameter_s22(r,theta,s22,s33,Xtildei,Utildei,Wtildei,n,beta);
updated_parameter_s22_never <- function(r,theta,s22,s33,qq,mmi,n){

# Do the Metropolis Step for s22.
# r = current value
# theta = current value
# s22 = current value
# s33 = current value

s22curr       <- s22
s22cand       <- s22 + (0.4 * (runif(1) - 0.5))
#GofSigmaecurr = formGofSigmae(r,theta,s22curr,s33,Xtildei,Utildei,Wtildei,n,beta);
#GofSigmaecand = formGofSigmae(r,theta,s22cand,s33,Xtildei,Utildei,Wtildei,n,beta);
# GofSigmaecurr = formGofSigmae(r,theta,s22curr,s33,qqa,qqb,qqc,qqd); #ab 3/5/11
# GofSigmaecand = formGofSigmae(r,theta,s22cand,s33,qqa,qqb,qqc,qqd);#ab 3/5/11

GofSigmaecurr = formGofSigmae_never_c(r,theta,s22curr,s33,qq,mmi)#ab 3/5/11
GofSigmaecand = formGofSigmae_never_c(r,theta,s22cand,s33,qq,mmi)#ab 3/5/11

gg            <- GofSigmaecand - GofSigmaecurr
gg            <- gg - ((mmi*n/2) * log(s22cand)) + ((mmi*n/2) * log(s22curr))

#ab 3/5/11
gg            <- min(1, exp(gg)*as.numeric((s22cand>=0)&(s22cand<=3)))
#gg            = min([1 exp(gg)]);
ss            <- pracma::rand(1,1)
s22new       <- (s22cand * as.numeric(ss < gg)) +  (s22curr * as.numeric(ss > gg))

return(s22new)

}

