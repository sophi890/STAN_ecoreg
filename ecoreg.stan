functions{ 
  real loglikeco_log(vector y, int numcounties, int[] numeffects, int[] numcats, matrix[] covlist, matrix adata, matrix whicha, vector pars){ //real
    vector[numeffects[1]] alphac = pars[1:(numeffects[1])]; // area-level covariate effects
    vector[sum(numcats)-numeffects[2]] alpha;
    vector[numeffects[3]] beta;
    vector[numcounties] N = adata[,2];
    vector[numcounties] q;
    real c = 16*sqrt(3) / (15 * pi());            
    vector[numcounties] denom;
    vector[prod(numcats)] cateffs = rep_vector(0.0, prod(numcats));
    //matrix[numcounties,prod(numcats)] alleffs;
    matrix[numcounties,prod(numcats)] plogis_alleffs;
    vector[prod(numcats)] ones = rep_vector(1.0, prod(numcats));
    vector[numcounties] p;
    vector[numcounties] loglik;
    
    // add group level covariate effect
    q = adata[,5:(numeffects[1]+4)] * alphac; // (3082 x 15)  x (15 x 1)
    
    if (numeffects[2]>0){ // if binary and categorical variables present
      alpha = pars[(numeffects[1]+1):(numeffects[1]+sum(numcats)-numeffects[2])]; // individual-level binary covariate effects
      cateffs = whicha*alpha; //cateffs by strata
    }
    else
      cateffs = rep_vector(0.0, prod(numcats));
    
    // add "normal" covariate effects
   if (numeffects[3]>0){
      beta = pars[(numeffects[1]+sum(numcats)-numeffects[2]+1):(numeffects[1]+sum(numcats)-numeffects[2]+numeffects[3])]; // individual-level normal covariate effects
      for (i in 1:num_elements(y)){
        denom[i] = sqrt(1+c*c*beta'*covlist[i]*beta);
      }
      q = q + (adata[,(numeffects[1]+5+prod(numcats)):(numeffects[1]+4+prod(numcats)+numeffects[3])] * beta) ./ denom; //adata[,404:405] numerator (3082 x 2) x (2 x 1). element wise division - both length 3082
    }
    
    
    // add categorical effects
    for (i in 1:num_elements(y)){
      for (j in 1:prod(numcats)){
        plogis_alleffs[i,j] = logistic_cdf(q[i] + cateffs[j],0,1);
      }
    }
   
    
    // end up with vector of probabilities for binomial, one for each FIPS 
    p = (adata[,(numeffects[1]+5):(numeffects[1]+4+prod(numcats))].* plogis_alleffs) * ones; //adata[,20:403]. (3082 x 384) x (384 x 1)
    
    // calculate loglikelihood
    for (i in 1:num_elements(y)){
      loglik[i] = y[i]*log(p[i]) + (N[i]-y[i])*log(1-p[i]);//binomial_lpmf(y[i] | N[i], p[i]); - not working b/c require ints not reals.
    }
    return sum(loglik);
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
} 

parameters{
  vector[30] pars; // use bounds? eg <lower=-3,upper=3>
}

model{
  pars[1] ~ normal(-6.5,5); // prior for intercept - based off of 0.0015=e^(-6.5) is approx baseline COVID death rate in US Feb
  pars[2:] ~ normal(0, sqrt(0.68)); // priors for remaining //corresponds to 95 percent odds ratio between 1/5 and 5
  y ~ loglikeco(numcounties, numeffects, numcats, covlist, adata, whicha, pars); // log likelihood
}

