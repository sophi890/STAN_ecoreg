library('rstan')

# Load data
load('adata.Rda')
load('whicha.Rda')
load('covlist.Rdata')

# Deaths
y = adata[,1]

# Run MCMC
fit = stan(file = 'ecoreg.stan', data = list(y = y, adata = adata, covlist = covlist, whicha = whicha, numeffects = c(15,13,2)))
# Getting log probabilitiy evalutes to log(0)! Error - need to debug something...

print(fit)
plot(fit)