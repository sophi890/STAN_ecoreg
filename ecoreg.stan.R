library('rstan')

# Load data
load('adata.Rda')
load('whicha.Rda')
load('covlist.Rdata')

# Deaths
y = adata[,1]

# Run MCMC
fit = stan(file = 'ecoreg.stan', data = list(y = y, adata = adata, covlist = covlist, whicha = whicha, numeffects = c(15,13,2)))

print(fit)

# Plot fit
plot(fit, show_density = TRUE, ci_level = 0.5, fill_color = "purple")
plot(fit, pars = c('pars[11]','pars[12]', 'pars[13]', 'pars[14]', 'pars[15]', 'pars[16]', 'pars[17]', 'pars[18]', 'pars[19]', 'pars[20]'), show_density = TRUE, ci_level = 0.5, fill_color = "purple")
plot(fit, pars = c('pars[21]','pars[22]', 'pars[23]', 'pars[24]', 'pars[25]', 'pars[26]', 'pars[27]', 'pars[28]', 'pars[29]', 'pars[30]'), show_density = TRUE, ci_level = 0.5, fill_color = "purple")

plot(fit, plotfun = 'trace', inc_warmup = TRUE)
traceplot(fit)

## Testing function - here loglik excludes binomial coefficients
expose_stan_functions(stanmodel = 'ecoreg.stan')
loglikeco_log(y=y, adata=adata, pars = rep(0, 30), covlist = covlist, whicha = whicha, numeffects = c(15,13,2))
