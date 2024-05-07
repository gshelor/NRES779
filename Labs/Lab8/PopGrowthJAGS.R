## PopGrowthJags.R

model{
  # priors
  K ~ dunif(0, 4000)
  r ~ dunif (0, 2)
  sigma ~ dunif(0, 2)
  tau <- 1/sigma^2
  # likelihood
  for(i in 1:n){
    mu[i] <- r*(1 - x[i]/K)
    y[i] ~ dnorm(mu[i], tau)
  }
  ## derived quantities for pop growth rate
  maxgrowth <- K/2
  for(j in 1:length(N)){
    growthrate[j] <- r * N[j] * (1 - N[j] / K)
  }
}