functions{ 
  real loglikeco_log(vector y, matrix adata, vector pars,  matrix[] covlist, matrix whicha, int[] numeffects){ //real
    vector[numeffects[1]] alphac = pars[1:(numeffects[1])]; // area-level covariate effects
    vector[numeffects[2]] alpha;
    vector[numeffects[3]] beta;
    vector[num_elements(y)] N = adata[,2];
    vector[num_elements(y)] q;
    real c = 16*sqrt(3) / (15 * pi());            
    vector[num_elements(y)] denom;
    vector[rows(whicha)] cateffs = rep_vector(0.0, rows(whicha));
    //matrix[num_elements(y),rows(whicha)] alleffs;
    matrix[num_elements(y),rows(whicha)] plogis_alleffs;
    vector[rows(whicha)] ones;
    vector[num_elements(y)] p;
    vector[num_elements(y)] loglik;
    
    // add group level covariate effect
    q = adata[,5:(numeffects[1]+4)] * alphac; // (3082 x 15)  x (15 x 1)
    
    // add "normal" covariate effects
   if (numeffects[3]>0){
      beta = pars[(numeffects[1]+numeffects[2]+1):(numeffects[1]+numeffects[2]+numeffects[3])]; // individual-level normal covariate effects
      for (i in 1:num_elements(y)){
        denom[i] = sqrt(1+c*c*beta'*covlist[i]*beta);
      }
      q = q + (adata[,(numeffects[1]+5+rows(whicha)):(numeffects[1]+4+rows(whicha)+numeffects[3])] * beta) ./ denom; //adata[,404:405] numerator (3082 x 2) x (2 x 1). element wise division - both length 3082
    }
    
    if (numeffects[2]>0){
      alpha = pars[(numeffects[1]+1):(numeffects[1]+numeffects[2])]; // individual-level binary covariate effects
      cateffs = whicha*alpha; //cateffs by strata
    }
    
    // add categorical effects
    for (i in 1:num_elements(y)){
      for (j in 1:rows(whicha)){
        plogis_alleffs[i,j] = logistic_cdf(q[i] + cateffs[j],0,1);
      }
    }
   
    ones = rep_vector(1.0, rows(whicha));
    
    // end up with vector of probabilities for binomial, one for each FIPS 
    p = (adata[,(numeffects[1]+5):(numeffects[1]+4+rows(whicha))].* plogis_alleffs) * ones; //adata[,20:403]. (3082 x 384) x (384 x 1)
    
    // calculate loglikelihood
    for (i in 1:num_elements(y)){
      loglik[i] = y[i]*log(p[i]) + (N[i]-y[i])*log(1-p[i]);//binomial_lpmf(y[i] | N[i], p[i]); - not working b/c require ints not reals.
    }
    return sum(loglik);
  }
}

data{ // how to make generalizable?
  vector[3082] y; // deaths
  matrix[3082,405] adata; // population + grouplevel covariates + normal covariates + percent ppl in each strata
  int numeffects[3];   // # group level, # categorical covariates, # normal covariates in that order
  matrix[numeffects[3],numeffects[3]] covlist[3082]; // covariance list for normal covariates
  matrix[384,13] whicha; // matrix of 1s and zeros based off of Jackson's whicha containing combinations (which categorical covariates a strata has)
} 

parameters{
  vector<lower=-3,upper=3>[30] pars; // use bounds?
}

model{
  //pars[1] ~ normal(-6.5,sqrt(5)); // prior for intercept - based off of 0.0015=e^(-6.5) is approx baseline COVID death rate in US
  pars ~ normal(0, sqrt(0.68)); // priors for remaining 
  y ~ loglikeco(adata, pars, covlist, whicha, numeffects); // log likelihood
}

