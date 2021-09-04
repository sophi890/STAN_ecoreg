functions{ 
  real loglikeco_log(vector y, int numcounties, int[] numeffects, matrix[] covlist, matrix adata, int[] states, vector cateffs, vector pars, vector rand){ //real
    vector[numeffects[1]] alphac = pars[1:(numeffects[1])]; // area-level covariate effects
    vector[numeffects[2]] beta; // normal covariate effects
    vector[numcounties] N = adata[,2]; // population each county
    vector[numcounties] q; 
    real c2 = 0.34584297349917959319510856775587; //  c = 16*sqrt(3) / (15 * pi()); 
    vector[numcounties] denom = rep_vector(1.0, numcounties); 
    vector[numcounties] normvec = rep_vector(0.0, numcounties); 
    matrix[numcounties, num_elements(cateffs)] plogis_alleffs;
    vector[num_elements(cateffs)] ones = rep_vector(1.0, num_elements(cateffs));
    vector[numcounties] p;

    // add group level covariate effect
    q = adata[,5:(numeffects[1]+4)] * alphac; // (3082 x 15)  x (15 x 1)

    
    // add "normal" covariate effects
   if (numeffects[2]>0){
      beta = pars[(numeffects[1]+1):(numeffects[1]+numeffects[2])]; // individual-level normal covariate effects
      for (i in 1:numcounties){
        denom[i] = sqrt(1+c2*beta'*covlist[i]*beta); // incorporates covariance 
      }
    normvec = adata[,(numeffects[1]+5+num_elements(cateffs)):(numeffects[1]+4+num_elements(cateffs)+numeffects[2])] * beta; 
    }

  for (j in 1:num_elements(cateffs)){ // progress along cols
    for (i in 1:numcounties){ // progress along rows within each col
      //plogis_alleffs[i,j] = 1 / (1 + exp(-1.0 * ((q[i] + cateffs[j] + normvec[i] + rand[states[i]]) ./ denom[i] )));
      plogis_alleffs[i,j] = logistic_cdf((q[i] + cateffs[j] + normvec[i] + rand[states[i]])/denom[i],0,1); 
    }
  }
    // end up with vector of probabilities for binomial, one for each FIPS 
    p = (adata[,(numeffects[1]+5):(numeffects[1]+4+num_elements(cateffs))].* plogis_alleffs) * ones; //adata[,20:403]. (3082 x 384) x (384 x 1)
    
    // calculate loglikelihood
    return sum(y .* log(p) + (N-y) .* log(1-p));
  }
}

data{ // how to make generalizable?
  int numcounties;
  vector[numcounties] y; // deaths
  vector[384] cateffs;
  int numeffects[2];   // # group level, # categorical covariates, # normal covariates in that order
  matrix[numeffects[2],numeffects[2]] covlist[numcounties]; // covariance list for normal covariates
  matrix[numcounties,numeffects[1]+numeffects[2]+num_elements(cateffs)+4] adata; // population + grouplevel covariates + normal covariates + percent ppl in each strata
  int states[numcounties];
} 

parameters{
  vector[numeffects[1]+numeffects[2]] pars; // use bounds? eg <lower=-3,upper=3>
  real<lower=0,upper=10> sigma;
  vector[49] rand;
}

model{
  rand ~ normal(0, sigma);
  pars[1] ~ normal(-10,5); // prior for intercept - based off of 0.0015=e^(-6.5) is approx baseline COVID death rate in US Feb
  pars[2:] ~ normal(0, sqrt(0.68)); // priors for remaining //corresponds to 95 percent odds ratio between 1/5 and 5
  y ~ loglikeco(numcounties, numeffects, covlist, adata, states, cateffs, pars, rand); // log likelihood
}

