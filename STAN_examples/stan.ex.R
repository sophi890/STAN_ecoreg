library("rstan") 
schools_dat <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

fit <- stan(file = 'STAN_examples/schools.stan', data = schools_dat)

print(fit)
plot(fit)
pairs(fit, pars = c("mu", "tau", "lp__"))

la <- extract(fit, permuted = TRUE) # return a list of arrays 
mu <- la$mu 

### return an array of three dimensions: iterations, chains, parameters 
a <- extract(fit, permuted = FALSE) 

### use S3 functions on stanfit objects
a2 <- as.array(fit)
m <- as.matrix(fit)
d <- as.data.frame(fit)

############### SECOND EXAMPLE EXPONENTIAL ##################

lambda = 1
dat = rexp(500, lambda)
len = length(dat)
fit = stan(file = 'STAN_examples/exponential.ex.stan', data = list(Y=dat, LENGTH=len))

print(fit)
mcmc_chain = as.matrix(fit)
acf(mcmc_chain[,'lambda'])
traceplot(fit)
pred = mcmc_chain[,'pred']
hist(dat, probability = T)
lines(density(pred))

plot(density(mcmc_chain[,'lambda']))

