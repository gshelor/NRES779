//
// This Stan program defines a model with a
// vector of values 'y' modeled as lognormally distributed
// with mean 'mu' and standard deviation 'sigma'.
//


// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] y; // SWE values
  vector[N] x1; // elevation
  vector[N] x2; // cumulative precipitation
  vector[N] x3; // temp min
  vector[N] x4; // temp mean
}


// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real b0; // intercept
  real b1; // elevation
  real b2; // cumulative precip
  real b3; // temp min
  real b4; // temp mean
  // real mu; // product of deterministic function, used in moment matching to determine shape and rate
  // real <lower=0> shape; // parameter in gamma dist
  // real <lower=0> rate; // parameter in gamma dist
  real<lower=0> sigma;
}

transformed parameters {
  real<lower=0> mu [N];
  real<lower=0> shape [N];
  real<lower=0> rate [N];
  for (i in 1:N){
    mu[i] = exp(b0 + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x4[i]);
    shape[i] = mu[i]^2 / sigma^2;
    rate[i] = mu[i] / sigma^2;
  }
}


// The model to be estimated.
// gamma dists using the shape and rate, not scale
model {
  // priors
  sigma ~ gamma(0.1,0.1);
  b0 ~ normal(0,10);
  b1 ~ normal(0,10);
  b2 ~ normal(0,10);
  b3 ~ normal(0,10);
  b4 ~ normal(0,10);
  // likelihood
  y ~ gamma(shape, rate);
}



generated quantities {
  int<lower=0, upper=1> mean_gt;
  int<lower=0, upper=1> sd_gt;
  vector[N] log_lik;
  {
    array[N] real y_rep = gamma_rng(shape, rate);
    mean_gt = mean(y_rep) > mean(y);
    sd_gt = sd(y_rep) > sd(y);
    for (i in 1:N) log_lik[i] = gamma_lpdf(y[i] | shape[i], rate[i]);
  }
}
