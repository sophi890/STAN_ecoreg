library('rstan')
options(mc.cores = 48)
rstan_options(auto_write = TRUE)

load("../data/ecoreg_random_96strata.RData")

fit96.random = stan(file = '../ecoreg_random.stan', 
                   data = list(y=adata[,1], 
                               numcounties = 3082, 
                               numeffects = c(15, 6, 2), 
                               numcats = c(2,2,2,2,2,3), 
                               covlist=covlist, 
                               adata=adata, 
                               whicha=whicha96, 
                               states=states), iter = 2000)


save(fit96.random,file=paste0("../fit96.random.Rdata"))
