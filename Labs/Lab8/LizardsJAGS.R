## LizardsJags.R

model{
  # priors
  b0 ~ dnorm(0, 1/4)
  b1 ~ dnorm(0, 1/4)
  # likelihood
  for(i in 1:n){
    logit(p[i]) <- b0 + b1*x[i]
    y[i] ~ dbern(p[i])
  }
  # derived quantities
  for(j in 1:60){
    logit(prob[j]) <- b0 + b1*OneToSixty_scale[j]
  }
}