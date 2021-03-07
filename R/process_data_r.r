 process_data <- function(rec_food, rec_energy, n, mmi, lambda_rec_food,
                          lambda_rec_energy){
###############################################################################
# This program process the data ready for analyses.
# INPUTS:
#     zz:                   Data set.
#     n:                    Sampel size.
#     mmi:                  Number of recalls want to use.
#     rec_food:             The recalls for the food
#     rec_energy:           The recalls for energy
#     lambda_rec_food:      Box-Cox transformation parameter for the food.
#     lambda_rec_energy:    Box-Cox transformation parameter for energy.
# OUTPUTS:
#     Wistar:               Indictor recall food > 0, n by mmi.
#     Wi2:                  Transformd and standardized the recalls foods, n by mmi.
#     Wi3:                  Transformed and standardized recall energy, n by mmi.
#     nointake_ffq:         Indicator people report no intake on ffq, n by 1.
#     didconsume:           Indicator people consume food on any of the recalls, n by 1.
#     a0_food:              Half of the minimum positive food value.
#     a0_energy:            Half of minimum energy value.
#     mumu, sigsig:         Mean and standard deviation of positive recall food values.
#     mu_e, sig_e:          Mean and standard deviation of recall energy.
###########################################################################
###########################################################################
# Get half of the minimum positive food value (a0_food) and half of minimum
# energy value (a0_energy)
###########################################################################
# library(R.matlab)
# library(pracma)
qq <- matrix(rec_food, nrow = n*mmi, ncol = 1, byrow = F)
uu <- as.numeric(qq > 0)
a0_food <- min(qq[uu == 1, ]) / 2
a0_energy <- min(matrix(rec_energy, nrow = n*mmi, ncol = 1, byrow = F)) /2
###########################################################################
# Transform and standardize recall energy
###########################################################################
Wi3           <- boxcoxtrans1(rec_energy, lambda_rec_energy)
mu_e           <- mean(matrix(Wi3, nrow = n*mmi, ncol = 1, byrow = F))
sig_e          <- sd(matrix(Wi3, nrow = n*mmi, ncol = 1, byrow = F))
Wi3           <- sqrt(2) * (Wi3 - mu_e) / sig_e
###########################################################################
# Transform and standardize the recalls for the foods.
# Make sure zeros stay as zeros!
###########################################################################
Wistar        <- matrix(as.numeric(rec_food > 0), n, mmi, byrow = F)
Wi2           <- boxcoxtrans1(rec_food, lambda_rec_food)
temp          <- matrix(Wi2, nrow = n*mmi, ncol = 1, byrow = F)
tempstar      <- matrix(Wistar, nrow = n*mmi, ncol = 1, byrow = F)
temp1         <- temp[tempstar == 1]
mumu          <- mean(temp1)
sigsig        <- sd(temp1)
Wi2           <- Wistar * sqrt(2) * (Wi2 - mumu) / sigsig
# indicator people consume food on any of the recalls
didconsume <- ifelse(rowSums(Wistar)>0, 1, 0)

return(list(Wistar = Wistar, Wi2 = Wi2, Wi3 = Wi3, didconsume = didconsume,
       a0_food = a0_food, a0_energy = a0_energy, mumu = mumu,
       sigsig = sigsig, mu_e = mu_e, sig_e = sig_e))
}
