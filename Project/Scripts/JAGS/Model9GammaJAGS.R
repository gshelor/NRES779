## modeling SWE with gamma dist in jags

model{
  # priors
  b0 ~ dnorm(0, 1/5)
  b1 ~ dnorm(0, 1/5)
  b2 ~ dnorm(0, 1/5)
  b3 ~ dnorm(0, 1/5)
  b4 ~ dnorm(0, 1/5)
  sigma ~ dgamma(1, 1/10)
  # likelihood
  for(i in 1:n){
    mu[i] <- b0 + b1*x1[i] + b2*x2[i] + b3*x3[i] + b4*x4[i]
    shape[i] <- mu[i]^2 / sigma^2
    rate[i] <- mu[i] / sigma^2
    y[i] ~ dgamma(shape[i], rate[i])
  }
}