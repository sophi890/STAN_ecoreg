library('rstan')

# Load data
load('data/adata.Rda')
load('data/whicha.Rda')
load('data/covlist.Rdata')

# Deaths
y = adata[,1]

# Run MCMC
fit1 = stan(file = 'ecoreg.stan', data = list(y = y, adata = adata, covlist = covlist, whicha = whicha, numeffects = c(15,13,2)), iter = 500)

print(fit1)

# Plot fit
plot(fit1, show_density = TRUE, ci_level = 0.5, fill_color = "purple")
plot(fit1, pars = c('pars[11]','pars[12]', 'pars[13]', 'pars[14]', 'pars[15]', 'pars[16]', 'pars[17]', 'pars[18]', 'pars[19]', 'pars[20]'), show_density = TRUE, ci_level = 0.5, fill_color = "purple")
plot(fit1, pars = c('pars[21]','pars[22]', 'pars[23]', 'pars[24]', 'pars[25]', 'pars[26]', 'pars[27]', 'pars[28]', 'pars[29]', 'pars[30]'), show_density = TRUE, ci_level = 0.5, fill_color = "purple")

plot(fit1, plotfun = 'trace', inc_warmup = TRUE)
traceplot(fit1)

## Testing function - here loglik excludes binomial coefficients - yes this matches Jackson's!
expose_stan_functions(stanmodel = 'ecoreg.stan')
loglikeco_log(y=y, adata=adata, pars = rep(0, 30), covlist = covlist, whicha = whicha, numeffects = c(15,13,2)) #-218254980 
loglikeco_log(y=y, adata=adata, pars = rep(-0.001, 30), covlist = covlist, whicha = whicha, numeffects = c(15,13,2)) #-203338046

# Jackson results: (very close - difference of C++ and R numerics?)
# -217301361
#-202384427