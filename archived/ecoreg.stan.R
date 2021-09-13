library('rstan')
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Model 1: Prelim example (simple ecological regression. Takes 3 hours to run 2000 iterations, converges.)

# Load data
load('data/adata_state_countyonly.Rda')

# Test loglik function
expose_stan_functions(stanmodel = 'countyonly.stan') # many warnings
loglikeco_log(y=adata[,1], numcounties=3082, numeffects=68, adata=adata, pars = rep(0, 68)) #-218254980
loglikeco_log(y=adata[,1], numcounties=3082, numeffects=68, adata=adata, pars = rep(-0.1, 68)) #-1269543
loglikeco_log(y=adata[,1], numcounties=3082, numeffects=68, adata=adata, pars = rep(0.1, 68)) #-2279137969

# Run MCMC
fit1 = stan(file = 'countyonly.stan', data = list(y=adata[,1], numcounties=3082, numeffects=68, adata=adata))
print(fit1)
plot(fit1, show_density = TRUE, ci_level = 0.5, fill_color = "purple")
plot(fit1, pars = c('pars[11]','pars[12]', 'pars[13]', 'pars[14]', 'pars[15]', 'pars[16]', 'pars[17]', 'pars[18]', 'pars[19]', 'pars[20]'), show_density = TRUE, ci_level = 0.5, fill_color = "purple")
plot(fit1, plotfun = 'trace', inc_warmup = TRUE)
plot(fit1, pars = c('pars[11]','pars[12]', 'pars[13]', 'pars[14]', 'pars[15]', 'pars[16]', 'pars[17]', 'pars[18]', 'pars[19]', 'pars[20]'), plotfun = 'trace', inc_warmup = TRUE)


# Model 2: Complex example 
#(includes joint distributions of categorical/binary covariates, normally distributed financial variables)

# Load data
load('data/adata.Rda')
load('data/whicha.Rda')
load('data/covlist.Rdata')

## Testing loglik function - yes this matches Jackson's! Although here I drop the binomial coeffs for simplicity.
expose_stan_functions(stanmodel = 'ecoreg.stan')
loglikeco_log(y=adata[,1], numcounties = 3082, numeffects = c(63, 6, 2), numcats = c(2,2,2,8,2,3), covlist=covlist, adata=adata, whicha=whicha, pars = rep(0, 78)) #-218254980  
loglikeco_log(y=adata[,1], numcounties = 3082, numeffects = c(63, 6, 2), numcats = c(2,2,2,8,2,3), covlist = covlist, adata=adata, whicha=whicha, pars = rep(-0.1, 78)) #-1242674
loglikeco_log(y=adata[,1], numcounties = 3082, numeffects = c(63, 6, 2), numcats = c(2,2,2,8,2,3), covlist = covlist, adata=adata, whicha=whicha, pars = rep(0.1, 78)) #-3077917034

# Run MCMC
fit2 = stan(file = 'ecoreg.stan', data = list(y=adata[,1], numcounties = 3082, numeffects = c(63, 6, 2), numcats = c(2,2,2,8,2,3), covlist=covlist, adata=adata, whicha=whicha))

print(fit2)
plot(fit2, show_density = TRUE, ci_level = 0.5, fill_color = "purple")
plot(fit2, pars = c('pars[11]','pars[12]', 'pars[13]', 'pars[14]', 'pars[15]', 'pars[16]', 'pars[17]', 'pars[18]', 'pars[19]', 'pars[20]'), show_density = TRUE, ci_level = 0.5, fill_color = "purple")
plot(fit2, pars = c('pars[21]','pars[22]', 'pars[23]', 'pars[24]', 'pars[25]', 'pars[26]', 'pars[27]', 'pars[28]', 'pars[29]', 'pars[30]'), show_density = TRUE, ci_level = 0.5, fill_color = "purple")
plot(fit2, plotfun = 'trace', inc_warmup = TRUE)
plot(fit2, pars = c('pars[11]','pars[12]', 'pars[13]', 'pars[14]', 'pars[15]', 'pars[16]', 'pars[17]', 'pars[18]', 'pars[19]', 'pars[20]'), plotfun = 'trace', inc_warmup = TRUE)
plot(fit2, pars = c('pars[21]','pars[22]', 'pars[23]', 'pars[24]', 'pars[25]', 'pars[26]', 'pars[27]', 'pars[28]', 'pars[29]', 'pars[30]'), plotfun = 'trace', inc_warmup = TRUE)

# Model 3: Complex example but now PM25 normal
#(includes joint distributions of categorical/binary covariates, normally distributed financial variables and Pm25)

# Load data
load('data/adata.pm25.Rda')
load('data/whicha.Rda')
load('data/covlist.pm25.Rda')

## Testing loglik function - yes this matches Jackson's! Although here I drop the binomial coeffs for simplicity.
expose_stan_functions(stanmodel = 'ecoreg.stan')
loglikeco_log(y=adata[,1], numcounties = 3082, numeffects = c(62, 6, 3), numcats = c(2,2,2,8,2,3), covlist = covlist.pm25, adata = adata, whicha = whicha, pars = rep(0, 78)) #-218254980  
loglikeco_log(y=adata[,1], numcounties = 3082, numeffects = c(62, 6, 3), numcats = c(2,2,2,8,2,3), covlist = covlist.pm25, adata = adata, whicha = whicha, pars = rep(-0.1, 78)) #-1241791
loglikeco_log(y=adata[,1], numcounties = 3082, numeffects = c(62, 6, 3), numcats = c(2,2,2,8,2,3), covlist = covlist.pm25, adata = adata, whicha = whicha, pars = rep(0.1, 78)) #-3076026276

# Run MCMC
fit3 = stan(file = 'ecoreg.stan', data = list(y=adata[,1], numcounties = 3082, numeffects = c(62, 6, 3), numcats = c(2,2,2,8,2,3), covlist = covlist.pm25, adata = adata, whicha = whicha), iter = 2000)

print(fit3)

# Try Model 1 with random state effects

# Load data
load('data/adata_random_countyonly.Rda')
load('data/states.Rda')

# Run MCMC
fit1.random = stan(file = 'countyonly_random.stan', data = list(y=adata[,1], numcounties=3082, numeffects=20, adata=adata, states=states))
print(fit1.random)

# Try Model 2 with random state effects
load('data/states.Rda')
load('data/adata_random.Rda')
load('data/whicha.Rda')
load('data/covlist.Rdata')

fit2.random = stan(file = 'ecoreg_random.stan', data = list(y=adata[,1], numcounties = 3082, numeffects = c(15, 6, 2), numcats = c(2,2,2,8,2,3), covlist=covlist, adata=adata, whicha=whicha, states=states))
print(fit2.random)

## Delete Later
for (i in 1:20){
  print(noquote(paste(signif(exp(mean(extract(fit1.random)$pars[,i])),4), ' & ','(', signif(exp(mean(extract(fit1.random)$pars[,i])-1.96*sd(extract(fit1.random)$pars[,i])),4), ', ', signif(exp(mean(extract(fit1.random)$pars[,i])+1.96*sd(extract(fit1.random)$pars[,i])),4), ')', sep = '')))
}
for (i in 1:30){
  print(noquote(paste(signif(exp(mean(extract(fit3)$pars[,i])),4), ' & ','(', signif(exp(mean(extract(fit3)$pars[,i])-1.96*sd(extract(fit3)$pars[,i])),4), ', ', signif(exp(mean(extract(fit3)$pars[,i])+1.96*sd(extract(fit3)$pars[,i])),4), ')', sep = '')))
}

strata <- adata[,grep("stratum", colnames(adata))] # extract strata from dataset
e <- rep(0, ncol(strata))
for (i in 1:ncol(strata)) { # loop over strata
  dat <- as.data.frame(adata[,"N"] * strata[,i]) # generate num people w/in strata
  mod <- glm(adata[,"y"] ~ ., family = poisson(), data = dat) # regress on number of deaths
  e[i] <- exp(coef(mod)[2]) - 1 # exponentiate coefficient and subtract 1 to get estimated risk parameter
}
summary(e) # median risk = 0.000526 (0.05% chance of death), max = 0.0246 (2.5% chance of death!), 24 NAs (no one in strata)
# use log(e/(1-e)) as gamma_s for each strata
e[is.na(e)] <- 0

fit3.random = stan(file = 'ecoreg_random_fixedoffset.stan', data = list(y=adata[,1], numcounties = 3082, numeffects = c(15, 2), covlist=covlist, adata=adata, states=states, cateffs = e), iter = 100)
