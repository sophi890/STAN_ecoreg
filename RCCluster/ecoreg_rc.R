library('rstan')
options(mc.cores = 24)
rstan_options(auto_write = TRUE)

#save(adata, covlist, whicha, file = "../data/ecoreg.RData")
#load("ecoreg.RData")
load("../data/ecoreg.RData")

y = adata[,1]

fit2 = stan(file = '../ecoreg.stan',
            data = list(y = y, 
                        numcounties = 3082, 
                        numeffects = c(63, 6, 2), 
                        numcats = c(2,2,2,8,2,3), 
                        covlist = covlist, 
                        adata = adata, 
                        whicha = whicha),
            iter = 2000)

save(fit2,file=paste0("../fit2.Rdata"))
