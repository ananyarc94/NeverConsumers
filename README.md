# NeverConsumers

The package is developed based on the paper : "Measurement Error Models With Zero Inflation and Multiple Sources of Zeros, With Applications to Never-Consumers in Nutrition ". In many problems, not limited to nutrition, there are variables measured with error that are not subject to occasional or excess zeros.There are practical cases, however, when some individuals will always report zero. This occurs, for example, with alcohol consumed from alcoholic beverages, reported on a daily basis, or with red meat, fish and deep orange vegetable consumption. In the former case, the majority of American adults are frequent or occasional consumers of alcoholic beverages,but a fraction never consume alcoholic beverages. Thus, in measurements made on a single occasion, e.g., one day, there are two sources of zeros: the occasional zeros caused by episodic consumption, and the hard zeros caused by being a never-consumer of alcoholic beverages. In most studies, it is not feasible to obtain more than a small number of repeated measurements, e.g. 2-4, on any individual. If all such measurements are zero it is not possible to distinguish,at the individual level, whether the person is an episodic consumer or a never consumer. This is what makes the problem so difficult.In many instances, nutritionists are interested in the amount of a dietary component consumed relative to the amount of calories consumed, thus controlling for energy consumption.

We propose a new model to accomplish five things: (a) adjust for measurement error; (b) allow for episodic consumers; (c) allow for never-consumers; (d) includes an additional continuous variable such as energy (caloric) intake; and (e) overcomes the numerical difficulties discussed above.


This package can be used to perform an entire Bayesian analysis using MCMC on the above mentioned model which is discussed in more details in the  paper and also in the attached vignette. We need to have installed the packages Pracma and R.matlab in order for this package to work.

Installation 
==================

To install, run the code below. You may need to install some other libraries since the package is not on Cran.

```{r}
library(devtools)
devtools::install_github("ananya.rc94/NeverConsumersR")
```

To install and build the vignette, run:

```{r}
library(devtools)
devtools::install_github("ananya.rc94/NeverConsumersR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
```

Example
==================
As an example, we have simulated a dataset to mimic the structure of the Eating at America’s Table Study (EATS) with 20 male subjects and 4 recalls per subject. The example dataset is attached with this package. The dataset has three responses - an indicator of whether food has been consumed on that particular day, amount of food consumed and amount of energy obtained. We also have 4 recalls for all subjects and some covariates. If there is any missing data, we drop those rows. For our simulated data we have four covariates age, BMI, DHQ_alcohol and DHQ_energy. We first load the dataset and specify all the MCMC parameters.

```{r}
zz   <-  read.csv("simulated_data.csv", header = T)# Name of the data set you will load, include the file path
nMCMC  <- 50000# Number of MCMC iterations. The more the better
nburn  <- 10000# Size of the burn-in
nthin  <- 50# thining
ndist  <- 200# average of the last ndist MCMC steps (after thinning)
                              # to get the cdf of usual intake
```

Next we have to specify to the function where the columns of interest are in the dataset. 


```{r}


X_cols          <- 2:5 # Columns where the covariates are
Food_col        <- 6 # Column where the episodic variable is
Energy_col      <- 7 # Column where energy (continuous) is
ID_col          <- 1 # Column for ID
FFQ_food_col    <- 4 # Column where the FFQ for the food is. 
FFQ_energy_col  <- 5 # Column where the FFQ for energy is. 
with_covariates_ind <- 3# What to include as covariates in the
                              # ever consumer model
                                # 0: a column of ones.
                                # 1: a column of ones, the FFQ, and
                                #    the indicator that the FFQ=0.
                                # 2: a column of ones and the FFQ.
                                # 3: a column of ones and the indicator
                                #    that the FFQ=0.
n_gen <- 30# Number of realizations of usual intake to
                                # generate. Must be a positive integer,
mmi             <- 4# Number of recalls, integer

```

Finally we have to specify the other parameters with no default values. 

```{r}

lambda_rec_food   <- 0.42
lambda_rec_food   <- 0.53


if (length(FFQ_food_col) > 0){
    lambda_FFQ_food   <- 0.0
}
if (length(FFQ_energy_col) > 0){
    lambda_FFQ_energy   <- 0.03
}

```
Then we make the function call as:

```{r}
nC_object = neverConsumers(zz, mmi, X_cols, Food_col, Energy_col, ID_col, FFQ_food_col,
               FFQ_energy_col, with_covariates_ind, n_gen, lambda_rec_food, 
               lambda_rec_energy, lambda_FFQ_food, lambda_FFQ_energy, nthin = nthin, 
               nMCMC = nMCMC, nburn = nburn, ndist = ndist)

```

This nC_object gives us the estimates of mean, standard deviation and a 95% confidence interval for all the parameters which are explained in more details in the function documentation. The function also generates the graphs of all the posterior distributions to make a complete analysis.
