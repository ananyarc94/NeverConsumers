#function rnew = updated_parameter_r(r,theta,s22,s33,Xtildei,Utildei,Wtildei,n,beta);
updated_parameter_r_never <- function(r,theta,s22,s33,qq,mmi,n){

# Do the Metropolis Step for r.
# r = current value
# theta = current value
# s22 = current value
# s33 = current value

rpossible <- seq(-0.99,0.99,length.out = 41)
spacing   <- rpossible[2] - rpossible[1]
if(abs(round(r,2)) <= 0.99 ){
    rcurr = round(r,2)
}else if(round(r,2) > 0.99){
    rcurr = 0.99
}else{
    rcurr = -0.99
}


if (rcurr == -0.99){
    i = 1
    ss    <- pracma::randn(1,1)
    rcand <- (rcurr * as.numeric(ss <= 0.33)) + ((rcurr + spacing)
            *as.numeric (ss > 0.33) *as.numeric(ss <= 0.66)) +
            + ((rcurr + (2 * spacing)) *as.numeric(ss > 0.66))
    
   
}
if (rcurr == 0.99){
    i = 2
    ss    <- pracma::randn(1,1)
    rcand <- (rcurr * as.numeric(ss <= 0.33)) + ((rcurr - spacing)
            * as.numeric(ss > 0.33) * as.numeric(ss <= 0.66)) +
            + ((rcurr - (2 * spacing)) * as.numeric(ss > 0.66))
}
if (rcurr > -0.99 && rcurr < 0.99){
    i = 3
    ss    <- pracma::randn(1,1); ss
    rcand <- (rcurr * as.numeric(ss <= 0.33)) + ((rcurr + spacing)
            * as.numeric(ss > 0.33) * as.numeric(ss <= 0.66)) +
                + ((rcurr - spacing) * as.numeric(ss > 0.66))
    
}

GofSigmaecurr = formGofSigmae_never(rcurr,theta,s22,s33,qq,mmi)
GofSigmaecand = formGofSigmae_never(rcand,theta,s22,s33,qq,mmi)

gg            <- GofSigmaecand - GofSigmaecurr
gg            <- gg - ((mmi*n/2) * log(1 - (rcand ^ 2))) + ((mmi*n/2)
                                                  * log(1 - (rcurr
                                                  ^ 2)))
gg            <- min(1, exp(gg))
ss1            <- pracma::randn(1,1)
rnew          <- (rcand * (ss1 < gg)) + ((rcurr * (ss1 > gg)))
return(rnew)
}


