#function GofSigmae = formGofSigmae_never(r,theta,s22,s33,Xtildei,Utildei,Wtildei,n,beta);

formGofSigmae_never <- function(r,theta,s22,s33,qq,mmi){

# Compute the loglikelihood contribution for Sigmae other than its
# determinant


  R     <-  matrix(c(1, 0, r*cos(theta),
                     0,  1, r*sin(theta),
                     r*cos(theta), r*sin(theta), 1), nrow = 3, byrow = T)
  
  A <- diag(c(1, sqrt(s22), sqrt(s33)))
  Sigmae       <- A %*% R %*% A
  iSigmae      <- pracma::inv(Sigmae)
  

tempMat <- matrix(0,ncol(qq),ncol(qq))
for (i in 1:mmi){
    tempMat <- tempMat + crossprod(qq[ , ,i])
}
tempMat <- iSigmae*tempMat
GofSigmae <- -0.5*(sum(tempMat))

return(GofSigmae = GofSigmae)
}
