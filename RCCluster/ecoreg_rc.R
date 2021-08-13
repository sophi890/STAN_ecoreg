library('rstan')
options(mc.cores = 24)
rstan_options(auto_write = TRUE)

load("./data/ecoreg.RData")

y = adata[,1]

fit2 = stan(file = './ecoreg.stan',
            data = list(y = y, 
                        numcounties = 3082, 
                        numeffects = c(15, 6, 2), 
                        numcats = c(2,2,2,8,2,3), 
                        covlist = covlist, 
                        adata = adata, 
                        whicha = whicha),
            iter = 1000)

save(fit2,file=paste0("/n/home01/xwu1993/STAN_ecoreg/fit2.Rdata"))

#save(adata, covlist, whicha, file = "ecoreg.RData")
#load("ecoreg.RData")
