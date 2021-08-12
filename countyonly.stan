functions{ 
  real loglikeco_log(vector y, matrix adata, vector pars){ //real
    vector[20] alphac = pars[1:20]; // area-level covariate effects
    vector[3082] N = adata[,2];
    vector[3082] q;
    vector[3082] p;
    vector[3082] loglik;
    
    // add group level covariate effect
    q = adata[,5:(20+4)] * alphac; // (3082 x 15)  x (15 x 1)
  
    for (i in 1:num_elements(y)){
        p[i] = logistic_cdf(q[i],0,1);
    }
    
    // end up with vector of probabilities for binomial, one for each FIPS 
    // calculate loglikelihood
    for (i in 1:num_elements(y)){
      loglik[i] = y[i]*log(p[i]) + (N[i]-y[i])*log(1-p[i]);//binomial_lpmf(y[i] | N[i], p[i]); - not working b/c require ints not reals.
    }
    return sum(loglik);
  }
}

data{ // how to make generalizable?
  vector[3082] y; // deaths
  matrix[3082,25] adata; // population + grouplevel covariates + normal covariates + percent ppl in each strata
} 

parameters{
  vector[20] pars; // use bounds? eg <lower=-3,upper=3>
}

model{
  pars[1] ~ normal(-6.5,5); // prior for intercept - based off of 0.0015=e^(-6.5) is approx baseline COVID death rate in US Feb
  pars[2:] ~ normal(0, sqrt(0.68)); // priors for remaining 
  y ~ loglikeco(adata, pars); // log likelihood
}

