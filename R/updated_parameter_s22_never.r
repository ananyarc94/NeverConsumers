#function s22new = updated_parameter_s22(r,theta,s22,s33,Xtildei,Utildei,Wtildei,n,beta);
updated_parameter_s22_never <- function(r_1,theta_1,s22_1,s33_1,qq_1,mmi_1,n_1){

# Do the Metropolis Step for s22.
# r_1 = current value
# theta_1 = current value
# s22_1 = current value
# s33_1 = current value

s22curr       <- s22_1
s22cand       <- s22_1 + (0.4 * (pracma::rand(1,1) - 0.5))
print(paste("this is s22curr = ", s22curr, " and s22cand = ", s22cand))
#GofSigmaecurr = formGofSigmae(r,theta,s22curr,s33,Xtildei,Utildei,Wtildei,n,beta);
#GofSigmaecand = formGofSigmae(r,theta,s22cand,s33,Xtildei,Utildei,Wtildei,n,beta);
# GofSigmaecurr = formGofSigmae(r,theta,s22curr,s33,qqa,qqb,qqc,qqd); #ab 3/5/11
# GofSigmaecand = formGofSigmae(r,theta,s22cand,s33,qqa,qqb,qqc,qqd);#ab 3/5/11

GofSigmaecurr = formGofSigmae_never(r_1,theta_1,s22curr,s33_1,qq_1,mmi_1)#ab 3/5/11
GofSigmaecand = formGofSigmae_never(r_1,theta_1,s22cand,s33_1,qq_1,mmi_1)#ab 3/5/11

gg            <- GofSigmaecand - GofSigmaecurr
gg            <- gg - ((mmi_1*n_1/2) * log(s22cand)) + ((mmi_1*n_1/2) * log(s22curr))
print(paste("gg initial", gg))
#ab 3/5/11
gg            <- min(1, exp(gg)*as.numeric((s22cand>=0)&(s22cand<=3)))
#gg            = min([1 exp(gg)]);
ss            <- pracma::rand(1,1)
print(paste("gg final = ", gg, " and ss = ", ss))
s22new       <- (s22cand * as.numeric(ss < gg)) +  (s22curr * as.numeric(ss > gg))

return(s22new)

}

# s22_r = vector()
# s22_c = vector()
# 
# for(i in 1:10000){
#   s22_r[i] = updated_parameter_s22_never(r,theta,s22,s33,qq,mmi,n)
#   s22_c[i] = updated_parameter_s22_never_c(r,theta,s22,s33,qq,mmi,n)
# }
# 
# mean(s22_r)
# mean(s22_c)
# 
# data.frame(s22_r, s22_c)
