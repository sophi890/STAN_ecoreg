library('rstan')
options(mc.cores = 24)
rstan_options(auto_write = TRUE)

load("../data/ecoreg.pm25.RData")

y = adata[,1]

fit3 = stan(file = '../ecoreg.stan',
            data = list(y = y, 
                        numcounties = 3082, 
                        numeffects = c(62, 6, 3), 
                        numcats = c(2,2,2,8,2,3), 
                        covlist = covlist.pm25, 
                        adata = adata, 
                        whicha = whicha),
            iter = 2000)

save(fit3,file=paste0("../fit3.Rdata"))
