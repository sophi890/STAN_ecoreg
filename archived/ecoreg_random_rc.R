library('rstan')
options(mc.cores = 24)
rstan_options(auto_write = TRUE)

load("../data/ecoreg_random.RData")

fit2.random = stan(file = '../ecoreg_random.stan', 
                   data = list(y=adata[,1], 
                               numcounties = 3082, 
                               numeffects = c(15, 6, 2), 
                               numcats = c(2,2,2,8,2,3), 
                               covlist=covlist, 
                               adata=adata, 
                               whicha=whicha, 
                               states=states))


save(fit2.random,file=paste0("../fit2.random.Rdata"))