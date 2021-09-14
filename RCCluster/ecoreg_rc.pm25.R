library('rstan')
options(mc.cores = 48)
rstan_options(auto_write = TRUE)

load("../data/ecoreg.pm25.RData")

fit3 = stan(file = '../ecoreg_random.stan',
            data = list(y = adata[,1], 
                        numcounties = 3082, 
                        numeffects = c(14, 6, 3), 
                        numcats = c(2,2,2,8,2,3), 
                        covlist = covlist.pm25, 
                        adata = adata, 
                        whicha = whicha,
                        states = states),
            iter = 2000)

save(fit3,file=paste0("../fit3.Rdata"))
