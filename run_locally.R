library('rstan')
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Model 1: Prelim example (simple ecological regression. Takes short time to run 2000 iterations, converges.)

## Load data
load('data/adata_random_countyonly.Rda')
load('data/states.Rda')

## Run MCMC
fit1.random = stan(file = 'countyonly_random.stan', data = list(y=adata[,1], numcounties=3082, numeffects=20, adata=adata, states=states))
print(fit1.random)
save(fit1.random, file = 'fit1.random.Rdata')



# Try Model 1 with random state effects exactly with Xiao's predictors

## Load data
load('data/adata_xiao.Rda')
load('data/states_x.Rda') # same states order

## Run MCMC
fit1_x.random = stan(file = 'countyonly_random.stan', data = list(y=adata[,1], numcounties=3089, numeffects=25, adata=adata, states=states))
print(fit1_x.random)
#save(fit1_x.random, file = 'fit1_x.random.RData')


# Model 2: Hierarchical
load('data/states.Rda')
load('data/adata_random.Rda')
load('data/whicha.Rda')
load('data/covlist.Rdata')
# save(adata, whicha, covlist, states, file = 'ecoreg_random.RData')

## This will not converge if you run locally (takes ages)
fit2.random = stan(file = 'ecoreg_random.stan', data = list(y=adata[,1], numcounties = 3082, numeffects = c(15, 6, 2), numcats = c(2,2,2,8,2,3), covlist=covlist, adata=adata, whicha=whicha, states=states))


# Model 3: Fixed offsets model

## Evaluate strata-specific risks

#strata <- adata[,grep("stratum", colnames(adata))] # extract strata from dataset
#e <- rep(0, ncol(strata))
#for (i in 1:ncol(strata)) { # loop over strata
  #dat <- as.data.frame(adata[,"N"] * strata[,i]) # generate num people w/in strata
  #mod <- glm(adata[,"y"] ~ ., family = poisson(), data = dat) # regress on number of deaths
  #e[i] <- exp(coef(mod)[2]) - 1 # exponentiate coefficient and subtract 1 to get estimated risk parameter
#}
#summary(e) # median risk = 0.000526 (0.05% chance of death), max = 0.0246 (2.5% chance of death!), 24 NAs (no one in strata)
# use log(e/(1-e)) as gamma_s for each strata
#gamma_s = log(e/(1-e)) # mean -7.133298
#gamma_s[is.na(gamma_s)] <- 0
# save(adata, covlist, whicha, states, gamma_s, file = 'ecoreg_fixedoffsets.RData')

## Load Data
load('data/ecoreg_fixedoffsets.RData')

## Run MCMC
fit3.random = stan(file = 'ecoreg_random_fixedoffset.stan', data = list(y=adata[,1], numcounties = 3082, numstrata = 384, numeffects = c(15, 2), covlist=covlist, adata=adata, states=states, cateffs = gamma_s))


# Model 3, condensing into 96 strata. Here age is collapsed from 8 categories into 2. 

## Load data
load('data/ecoreg_fixedoffsets_96strata.RData')

## Run MCMC
fit3.random.age2 = stan(file = 'ecoreg_random_fixedoffset.stan', data = list(y = adata[,1], numcounties = 3082, numstrata = 96, numeffects = c(15, 2), covlist = covlist, adata = adata, states = states, cateffs = gamma_s))
