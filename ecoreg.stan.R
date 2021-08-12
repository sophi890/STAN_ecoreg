library('rstan')
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Model 1: Prelim example (simple ecological regression. Takes 35 minutes to run 2000 iterations, converges.)

# Load data
load('data/adata_countyonly.Rda') 
y = adata[,1]

# Test loglik function
expose_stan_functions(stanmodel = 'countyonly.stan') # many warnings
loglikeco_log(y=y, adata=adata, pars = rep(0, 20)) #-218254980
loglikeco_log(y=y, adata=adata, pars = rep(-0.1, 20)) #-1293860

fit1 = stan(file = 'countyonly.stan', data = list(y=y, numcounties=3082, numeffects=20, adata=adata))
print(fit1)
plot(fit1, show_density = TRUE, ci_level = 0.5, fill_color = "purple")
plot(fit1, pars = c('pars[11]','pars[12]', 'pars[13]', 'pars[14]', 'pars[15]', 'pars[16]', 'pars[17]', 'pars[18]', 'pars[19]', 'pars[20]'), show_density = TRUE, ci_level = 0.5, fill_color = "purple")
plot(fit1, plotfun = 'trace', inc_warmup = TRUE)


# Model 2: Complex example 
#(includes joint distributions of categorical/binary covariates, normally distributed financial variables)

# Load data
load('data/adata.Rda')
load('data/whicha.Rda')
load('data/covlist.Rdata')

# Deaths
y = adata[,1]

## Testing loglik function - yes this matches Jackson's! Although here I drop the binomial coeffs for simplicity.
expose_stan_functions(stanmodel = 'ecoreg.stan')
loglikeco_log(y = y, numcounties = 3082, numeffects = c(15, 6, 2), numcats = c(2,2,2,8,2,3), covlist=covlist, adata=adata, whicha=whicha, pars = rep(0, 30)) #-218254980  
loglikeco_log(y=y, numcounties = 3082, numeffects = c(15, 6, 2), numcats = c(2,2,2,8,2,3), covlist = covlist, adata=adata, whicha=whicha, pars = rep(-0.1, 30)) #-1233863
loglikeco_log(y=y, numcounties = 3082, numeffects = c(15, 6, 2), numcats = c(2,2,2,8,2,3), covlist = covlist, adata=adata, whicha=whicha, pars = rep(0.1, 30)) #-3047040336

# Run MCMC
fit2 = stan(file = 'ecoreg.stan', data = list(y = y, numcounties = 3082, numeffects = c(15, 6, 2), numcats = c(2,2,2,8,2,3), covlist=covlist, adata=adata, whicha=whicha))

print(fit2)
plot(fit2, show_density = TRUE, ci_level = 0.5, fill_color = "purple")
plot(fit2, pars = c('pars[11]','pars[12]', 'pars[13]', 'pars[14]', 'pars[15]', 'pars[16]', 'pars[17]', 'pars[18]', 'pars[19]', 'pars[20]'), show_density = TRUE, ci_level = 0.5, fill_color = "purple")
plot(fit2, pars = c('pars[21]','pars[22]', 'pars[23]', 'pars[24]', 'pars[25]', 'pars[26]', 'pars[27]', 'pars[28]', 'pars[29]', 'pars[30]'), show_density = TRUE, ci_level = 0.5, fill_color = "purple")
plot(fit2, plotfun = 'trace', inc_warmup = TRUE)

# Model 3: Normal only example
#(includes only normally distributed financial variables)

# Load data
#load('data/adata3.Rda')
#load('data/covlist.Rdata')

# Deaths
#y = adata[,1]

# Run MCMC
#fit5 = stan(file = 'ecoreg_normalonly.stan', data = list(y = y, covlist = covlist, adata=adata))

#print(fit5)
#plot(fit5, show_density = TRUE, ci_level = 0.5, fill_color = "purple")
#plot(fit5, pars = c('pars[11]','pars[12]', 'pars[13]', 'pars[14]', 'pars[15]', 'pars[16]', 'pars[17]', 'pars[18]', 'pars[19]', 'pars[20]'), show_density = TRUE, ci_level = 0.5, fill_color = "purple")
#plot(fit5, plotfun = 'trace', inc_warmup = TRUE)

# Model 4: Only age/sex/race joint dist (from 2010.) normal variables too.

# Load data
#load('data/adata_agesexraceonly.Rda')
#load('data/covlist.Rdata')
#load('data/whicha_agesexraceonly.Rda')

# Deaths
#y = adata[,1]

# Run MCMC
#fit6 = stan(file = 'ecoreg.stan', data = list(y = y, numcounties = 3082, numeffects = c(18, 3, 2), numcats = c(8,2,3), covlist=covlist, adata=adata, whicha=whicha.mat))

#print(fit6)
#plot(fit6, show_density = TRUE, ci_level = 0.5, fill_color = "purple")
#plot(fit6, pars = c('pars[11]','pars[12]', 'pars[13]', 'pars[14]', 'pars[15]', 'pars[16]', 'pars[17]', 'pars[18]', 'pars[19]', 'pars[20]'), show_density = TRUE, ci_level = 0.5, fill_color = "purple")
#plot(fit6, plotfun = 'trace', inc_warmup = TRUE)
