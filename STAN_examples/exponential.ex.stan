functions{
  real newexp_log(vector x, real lam){
    vector[num_elements(x)] prob;
    real lprob;
    for (i in 1:num_elements(x)){
      prob[i] = lam*exp(-lam*x[i]);
    }
    lprob = sum(log(prob));
    return lprob;
  }
}

data {
  int LENGTH;
  vector[LENGTH] Y;
}

parameters{
  real<lower=0> lambda;
}

model{
  real alpha;
  real beta;
  alpha = 1.0;
  beta = 1.0;
  lambda ~ gamma(alpha, beta);
  Y ~ newexp(lambda);
}

generated quantities{
  real pred;
  pred = exponential_rng(lambda);
}