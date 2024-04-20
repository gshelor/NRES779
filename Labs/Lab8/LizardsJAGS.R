## LizardsJags.R

model{
  # priors
  b0 ~ dnorm(0, 1/10000)
  b1 ~ dnorm(0, 1/10000)
  # likelihood
  for(i in 1:n){
    logit(p[i]) <- b0 + b1*x[i]
    y[i] ~ dbern(p[i])
  }
}