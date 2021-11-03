functions{ 
  real loglikeco_log(vector y, int numcounties, int[] numeffects, int[] numcats, matrix[] covlist, matrix adata, matrix whicha, int[] states, vector pars, vector rand){ //real
    vector[numeffects[1]] alphac = pars[1:(numeffects[1])]; // area-level covariate effects
    vector[sum(numcats)-numeffects[2]] alpha; // binary/categorical covariate effects
    vector[numeffects[3]] beta; // normal covariate effects
    vector[numcounties] N = adata[,2]; // population each county
    vector[numcounties] q; 
    real c2 = 0.34584297349917959319510856775587; //  c = 16*sqrt(3) / (15 * pi()); 
    vector[numcounties] denom = rep_vector(1.0, numcounties); 
    vector[numcounties] normvec = rep_vector(0.0, numcounties);
    vector[prod(numcats)] cateffs = rep_vector(0.0, prod(numcats));
    vector[numcounties] p = rep_vector(0.0, numcounties);
    
    // add group level covariate effect
    q = adata[,5:(numeffects[1]+4)] * alphac; // (3082 x 15)  x (15 x 1)
    
    if (numeffects[2]>0){ // if binary and categorical variables present
      alpha = pars[(numeffects[1]+1):(numeffects[1]+sum(numcats)-numeffects[2])]; // individual-level binary covariate effects
      cateffs = whicha*alpha; //cateffs by strata
    }
    
    // add "normal" covariate effects
   if (numeffects[3]>0){
      beta = pars[(numeffects[1]+sum(numcats)-numeffects[2]+1):(numeffects[1]+sum(numcats)-numeffects[2]+numeffects[3])]; // individual-level normal covariate effects
      for (i in 1:numcounties){
        denom[i] = sqrt(1+c2*beta'*covlist[i]*beta);
      }
    normvec = adata[,(numeffects[1]+5+prod(numcats)):(numeffects[1]+4+prod(numcats)+numeffects[3])] * beta; 
    }

  for (j in 1:prod(numcats)){ // progress along cols
    p += adata[,(numeffects[1]+4+j)] ./ (1.0 + exp(-(cateffs[j] + q + normvec + rand[states]) ./ denom));
  }
    // end up with vector of probabilities for binomial, one for each FIPS 
    
    // calculate loglikelihood
    return sum(y .* log(p) + (N-y) .* log(1-p));
  }
}

data{ // how to make generalizable?
  int numcounties;
  vector[numcounties] y; // deaths
  int numeffects[3];   // # group level, # categorical covariates, # normal covariates in that order
  int numcats[numeffects[2]]; // contains # of levels for each categorical/binary variable eg c(8,2,3) for age, race, sex
  matrix[numeffects[3],numeffects[3]] covlist[numcounties]; // covariance list for normal covariates
  matrix[numcounties,numeffects[1]+numeffects[3]+prod(numcats)+4] adata; // population + grouplevel covariates + normal covariates + percent ppl in each strata
  matrix[prod(numcats),sum(numcats)-numeffects[2]] whicha; // matrix of 1s and zeros based off of Jackson's 'whicha' containing combinations (which categorical covariates a strata has)
  int states[numcounties];
} 

parameters{
  vector[numeffects[1]+sum(numcats)-numeffects[2]+numeffects[3]] pars; // use bounds? eg <lower=-3,upper=3>
  real<lower=0,upper=1> sigma;
  vector[49] rand;
}

model{
  rand ~ normal(0, sigma);
  pars[1] ~ normal(-8,5); // prior for intercept - based off of 0.0015=e^(-6.5) is approx baseline COVID death rate in US Feb
  pars[2:] ~ normal(0, sqrt(0.68)); // priors for remaining //corresponds to 95 percent odds ratio between 1/5 and 5
  y ~ loglikeco(numcounties, numeffects, numcats, covlist, adata, whicha, states, pars, rand); // log likelihood
}
