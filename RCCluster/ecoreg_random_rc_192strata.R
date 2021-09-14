library('rstan')
options(mc.cores = 48)
rstan_options(auto_write = TRUE)

load("../data/ecoreg_192strata.RData")

fit192.random = stan(file = '../ecoreg_random.stan', 
                   data = list(y=adata[,1], 
                               numcounties = 3082, 
                               numeffects = c(15, 6, 2), 
                               numcats = c(2,2,2,4,2,3), 
                               covlist=covlist, 
                               adata=adata, 
                               whicha=whicha192, 
                               states=states), iter = 1000)


save(fit192.random,file=paste0("../fit192.random.Rdata"))
