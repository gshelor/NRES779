//
// This Stan program defines a model with a
// vector of values 'y' modeled as lognormally distributed
// with mean 'mu' and standard deviation 'sigma'. It also features an interactive effect between mean temp 
// and daily precipitation.


// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] y; // SWE values, plus 0.001
  vector[N] x1; // temp mean
  vector[N] x2; // daily precip
}

// The parameters accepted by the model. Our model
// accepts parameters b0, b1, and 'sigma'.
parameters {
  real b0; // intercept
  real b1; // temp mean/precip interaction
  real<lower=0> sigma;
}

transformed parameters {
  real mu [N];
  for (i in 1:N){
    mu[i] = b0 + b1*x1[i]*x2[i];
  }
}

// The model to be estimated. We model the output
// 'y' to be lognormally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  b0 ~ normal(0,10);
  b1 ~ normal(0,10);
  sigma ~ gamma(1,1)
  y ~ lognormal(mu, sigma);
}

generated quantities {
  int<lower=0, upper=1> meanlog_gt;
  int<lower=0, upper=1> sdlog_gt;
  vector[N] log_lik;
  {
    array[N] real y_rep = lognormal_rng(mu, sigma);
    meanlog_gt = mean(y_rep) > mean(y);
    sdlog_gt = sd(y_rep) > sd(y);
    for (i in 1:N) log_lik[i] = lognormal_lpdf(y[i] | mu[i], sigma);
  }
}