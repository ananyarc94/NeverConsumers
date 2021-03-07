#function thetanew = updated_parameter_theta(r,theta,s22,s33,Xtildei,Utildei,Wtildei,n,beta);
updated_parameter_theta_never <- function(r,theta,s22,s33,qq,mmi){

# Do the Metropolis Step for theta.
# r = current value
# theta = current value
# s22 = current value
# s33 = current value

thetapossible <- pi * seq(-0.99,0.99, length.out = 41)
thetamin      <- min(thetapossible, theta)
thetamax      <- max(thetapossible, theta)
spacing       <- thetapossible[2] - thetapossible[1]
thetacurr     <- theta
if (thetacurr <= thetamin){
    ss        <- pracma::randn(1,1)
    thetacand <- (thetacurr * as.numeric(ss <= 0.33)) + ((thetacurr + spacing)
                * as.numeric(ss > 0.33) * as.numeric(ss <= 0.66)) +
                + ((thetacurr + (2 * spacing)) * as.numeric(ss > 0.66))
}
if (thetacurr >= thetamax){
    ss        <- pracma::randn(1,1)
    thetacand <- (thetacurr * as.numeric(ss <= 0.33)) + ((thetacurr - spacing)
                * as.numeric(ss > 0.33) * as.numeric(ss <= 0.66)) +
                + ((thetacurr - (2 * spacing)) * as.numeric(ss > 0.66))
}
if (thetacurr  > thetamin){
    if (thetacurr < thetamax){
    ss        <- pracma::randn(1,1)
    thetacand <- (thetacurr * as.numeric(ss <= 0.33)) + ((thetacurr + spacing)
                * as.numeric(ss > 0.33) * as.numeric(ss <= 0.66)) +
                + ((thetacurr - spacing) * as.numeric(ss > 0.66))
    }
}
GofSigmaecurr = formGofSigmae_never(r,thetacurr,s22,s33,qq,mmi)
GofSigmaecand = formGofSigmae_never(r,thetacand,s22,s33,qq,mmi)

gg            <- min(1, exp(GofSigmaecand - GofSigmaecurr))
ss            <- pracma::randn(1,1)
thetanew      <- (thetacand * as.numeric(ss < gg)) + (thetacurr * as.numeric(ss > gg))
thetanew      <- ifelse(thetanew > 1.2, 1, thetanew)
return(thetanew)
}
