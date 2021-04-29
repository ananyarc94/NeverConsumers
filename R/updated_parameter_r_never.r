#function rnew = updated_parameter_r(r,theta,s22,s33,Xtildei,Utildei,Wtildei,n,beta);
updated_parameter_r_never <- function(r,theta,s22,s33,qq,mmi,n){

# Do the Metropolis Step for r.
# r = current value
# theta = current value
# s22 = current value
# s33 = current value

rpossible <- seq(-0.99,0.99,length.out = 41)
spacing   <- rpossible[2] - rpossible[1]
rcurr = as.numeric(r)
# if(abs(round(r,2)) <= 0.99 ){
#     rcurr = round(r,2)
# }else if(round(r,2) > 0.99){
#     rcurr = 0.99
# }else{
#     rcurr = -0.99
# }
# if(r <= -0.99){
#     rcurr = 0.99
# }else if( r >= 0.99){
#     rcurr = 0.99
# }else{
#    rcurr = round(r,2) 
# }
# print(rcurr)
# print(typeof(rcurr))
# print(rcurr > -0.99 && rcurr < 0.99)

if(isTRUE(all.equal(rcurr, -0.99))){
    i = 1
    ss    <- pracma::rand(1,1);ss
    rcand <- (rcurr * as.numeric(ss <= 0.33)) + ((rcurr + spacing)
            *as.numeric (ss > 0.33) *as.numeric(ss <= 0.66)) +
            + ((rcurr + (2 * spacing)) *as.numeric(ss > 0.66))
    #print(paste("this is ss:", ss,"when rcurr = ",rcurr,"and rcand = ",rcand,"and spacing = ",spacing," and i = ", i))
    
   
} else if(isTRUE(all.equal(rcurr, 0.99))){
    i = 2
    ss    <- pracma::rand(1,1)
    rcand <- (rcurr * as.numeric(ss <= 0.33)) + ((rcurr - spacing)
            * as.numeric(ss > 0.33) * as.numeric(ss <= 0.66)) +
            + ((rcurr - (2 * spacing)) * as.numeric(ss > 0.66))
    #print(paste("this is ss:", ss,"when rcurr = ",rcurr,"and rcand = ",rcand,"and spacing = ",spacing," and i = ", i))
} else{
    # if(rcurr > -0.99){
    #     if(rcurr < 0.99){
    i = 3
    ss    <- pracma::rand(1,1); ss
    rcand <- (rcurr * as.numeric(ss <= 0.33)) + ((rcurr + spacing)
            * as.numeric(ss > 0.33) * as.numeric(ss <= 0.66)) +
                + ((rcurr - spacing) * as.numeric(ss > 0.66))
    #print(paste("this is ss:", ss,"when rcurr = ",rcurr,"and rcand = ",rcand,"and spacing = ",spacing," and i = ", i))
    
}

# if(rcand <= -0.99){
#     rcand = 0.99
# }else if( rcand >= 0.99){
#     rcand = 0.99
# }else{
#     rcand = round(rcand,2) 
# }

#print(paste("this is rcand =",rcand,"before we call formGOF"))
GofSigmaecurr = formGofSigmae_never(rcurr,theta,s22,s33,qq,mmi)
GofSigmaecand = formGofSigmae_never(rcand,theta,s22,s33,qq,mmi)

gg            <- GofSigmaecand - GofSigmaecurr
gg            <- gg - ((mmi*n/2) * log(1 - (rcand ^ 2))) + ((mmi*n/2)
                                                  * log(1 - (rcurr
                                                  ^ 2)))
gg = ifelse(gg < -150, -150, gg) 
gg            <- min(1, exp(gg))
ss1            <- pracma::rand(1,1)
rnew          <- (rcand * (ss1 < gg)) + ((rcurr * (ss1 > gg)))


#print(paste("this is ss1 = ", ss1, "and gg = ", gg," and rnew = ",rnew))

# if(abs(round(r,4)) <= 0.9999 ){
#     rnew = rnew
# }else if(round(r,4) > 0.9999){
#     rnew = 0.9999
# }else{
#     rnew = -0.9999
# }

return(rnew)
}

# g_r = g_c = vector()
# for(i in 1:1000) {
#   g_r[i] = formGofSigmae_never(rcurr,theta,s22,s33,qq,mmi)
#   g_c[i] = formGofSigmae_never_c(rcurr,theta,s22,s33,qq,mmi)
# }

# r_r = r_c = vector()
# for(i in 1:10000) {
#   r_r[i] = updated_parameter_theta_never(r,theta,s22,s33,qq,mmi)
#   r_c[i] = updated_parameter_theta_never_c(r,theta,s22,s33,qq,mmi)
# }

