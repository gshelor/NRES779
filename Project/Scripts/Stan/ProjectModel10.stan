//
// This Stan program defines a model with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] y; // SWE values, plus 0.001
  vector[N] x1; // elevation
  vector[N] x2; // cumulative precipitation
  vector[N] x3; // temp min
}

// The parameters accepted by the model. Our model
// accepts parameters b0,b1, b2, b3, and 'sigma'.
parameters {
  real b0; // intercept
  real b1; // elevation
  real b2; // cumulative precip
  real b3; // temp min
  real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'. I am also truncating the distribution with a lower limit of 0.
model {
  b0 ~ normal(0,10);
  b1 ~ normal(0,100);
  b2 ~ normal(0,10);
  b3 ~ normal(0,10);
  y ~ normal(b0 + b1*x1 + b2*x2 + b3*x3, sigma) T[0,];
}

generated quantities {
  int<lower=0, upper=1> mean_gt;
  int<lower=0, upper=1> sd_gt;
  vector[N] log_lik;
  {
    array[N] real y_rep = normal_rng(b0 + b1*x1 + b2*x2 + b3*x3, sigma);
    mean_gt = mean(y_rep) > mean(y);
    sd_gt = sd(y_rep) > sd(y);
    for (i in 1:N) log_lik[i] = normal_lpdf(y[i] | b0 + b1*x1[i] + b2*x2[i] + b3*x3[i], sigma);
  }
}
