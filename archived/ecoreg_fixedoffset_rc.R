library('rstan')
options(mc.cores = 48)
rstan_options(auto_write = TRUE)

load("../data/ecoreg_fixedoffsets.RData")

fit.fixedoffsets = stan(file = 'ecoreg_random_fixedoffset.stan', 
                   data = list(y=adata[,1], 
                               numcounties = 3082, 
                               numeffects = c(15, 2), 
                               covlist=covlist, 
                               adata=adata, 
                               states=states, 
                               cateffs = gamma_s),
                   iter = 1000)

save(fit.fixedoffsets,file=paste0("../fit_fixedoffsets.Rdata"))
