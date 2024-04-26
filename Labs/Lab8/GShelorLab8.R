## loading packages
library(pacman)
p_load(here, tidyverse, rjags, ESS575, LaplacesDemon, coda)


##### Section 2 Question 1 #####
## pulling dataset from ESS575 package
Logistic = Logistic[order(Logistic$PopulationSize),]

## creating vector of N values
N <- seq(0,1100,10)
## initializing model
inits = list(
  list(K = 1500, r = .2, sigma = 1),
  list(K = 1000, r = .15, sigma = .1),
  list(K = 900, r = .3, sigma = .01))
data_pop = list(
  n = nrow(Logistic),
  x = as.double(Logistic$PopulationSize),
  y = as.double(Logistic$GrowthRate),
  N = N)
n.adapt = 5000
n.update = 10000
n.iter = 10000
# Call to JAGS
set.seed(1)
## pop growth model
jm_pop = jags.model(here("Labs", "Lab8", "PopGrowthJAGS.R"), data = data_pop, inits = inits,
                    n.chains = length(inits), n.adapt = n.adapt)
update(jm_pop, n.iter = n.update)
zm_pop = coda.samples(jm_pop, variable.names = c("K", "r", "sigma", "tau", "mu"),
                                  n.iter = n.iter, n.thin = 1)
plot_pop = coda.samples(jm_pop, variable.names = c("K", "r", "sigma", "tau"),
                      n.iter = n.iter, n.thin = 1)
zj_pop = jags.samples(jm_pop, variable.names = c("K", "r", "sigma", "tau", "mu", "growthrate"), n.iter = n.iter, n.thin = 1)
## estimating posterior dist
post_df_pop = as.data.frame(rbind(zm_pop[[1]], zm_pop[[2]], zm_pop[[3]]))

## plotting densities of params
par(mfrow = c(2,2))
hist(post_df_pop$K, prob = TRUE)
hist(post_df_pop$r, prob = TRUE)
hist(post_df_pop$sigma, prob = TRUE)
hist(post_df_pop$tau, prob = TRUE)
## plotting trace plots and density plots
par(mfrow = c(2,4), mar = c(2, 2, 2, 2))
plot(plot_pop)


##### Section 2 Question 2 #####
## plotting median growth rate with credible intervals as a 
growthrate_summary = summary(zj_pop$growthrate, quantile, c(.025, .5, .975))$stat
par(mfrow = c(1,1), mar = c(4,4,4,4))
plot(N, growthrate_summary[2,], type = 'l', ylim = c(min(growthrate_summary[1,]), max(growthrate_summary[3,])), xlab = "N (Population)", ylab = "Growth Rate", main = "Median Growth Rate as a Function of Population, with 95% Credible Intervals")
lines(N, growthrate_summary[1,], lty= 2)
lines(N, growthrate_summary[3,], lty = 2)


##### Section 2 Question 3 #####
1 - ecdf(zj_pop$r)(0.22)
ecdf(zj_pop$r)(0.22) - ecdf(zj_pop$r)(0.18)


##### Section 3 Question 4 #####
IslandLizards <- IslandLizards
scaledPAR <- scale(IslandLizards$perimeterAreaRatio)
IslandLizards$scaledPAR <- scaledPAR[,1]


##### Section 3 Question 6 #####
## creating a vector of 1:60 for question 8
OneToSixty <- 1:60
## scaling it just as the PAR were in the collected data
OneToSixty_scale <- scale(OneToSixty)[,1]
## initializing model
inits_liz = list(
  list(b0 = 0, b1 = 2),
  list(b0 = 2, b1 = 5),
  list(b0 = 5, b1 = 0))
data_liz = list(
  n = nrow(IslandLizards),
  x = as.double(IslandLizards$scaledPAR),
  y = as.double(IslandLizards$presence),
  OneToSixty_scale = OneToSixty_scale)
n.adapt = 5000
n.update = 10000
n.iter = 10000
# Call to JAGS
set.seed(1)
## pop growth model
jm_liz = jags.model(here("Labs", "Lab8", "LizardsJAGS.R"), data = data_liz, inits = inits_liz,
                    n.chains = length(inits_liz), n.adapt = n.adapt)
update(jm_liz, n.iter = n.update)
zm_liz = coda.samples(jm_liz, variable.names = c("b0", "b1", "prob"),
                      n.iter = n.iter, n.thin = 1)
zj_liz = jags.samples(jm_liz, variable.names = c("b0", "b1", "prob"), n.iter = n.iter, n.thin = 1)
## converting coda object to df
post_df_liz = as.data.frame(rbind(zm_liz[[1]], zm_liz[[2]], zm_liz[[3]]))


##### Section 3 Question 7 #####
## plotting densities of params
par(mfrow = c(2,1))
hist(post_df_liz$b0, prob = TRUE)
hist(post_df_liz$b1, prob = TRUE)
## plotting trace plots and density plots
par(mfrow = c(2,3), mar = c(2, 2, 2, 2))
plot(zm_liz[,1:2])

## Gelman diagnostic
gelman.diag(zm_liz[,1:2])
heidel.diag(zm_liz[,1:2])


##### Section 3 Question 8 #####
prob_med <- apply(post_df_liz[,3:62], 2, median) #median
prob_lower <- apply(post_df_liz[,3:62], 2, quantile, 0.025) #lower CI
prob_upper <- apply(post_df_liz[,3:62], 2, quantile, 0.975) 
## plotting median growth rate with credible intervals as a dashed line
par(mfrow = c(1,1), mar = c(4,4,4,4))
plot(IslandLizards$perimeterAreaRatio, IslandLizards$presence, xlab = "Perimeter to Area Ratio", ylab = "Probability of Presence", main = "Probability of Lizard Presence as Function of PAR", sub = "Median of Predicted Probability = solid line, 95% credible intervals = dashed line")
lines(prob_med, lty = 1)
lines(prob_lower, lty = 2)
lines(prob_upper, lty = 2)

##### Section 3 Question 9 #####
diff <- post_df_liz$`prob[10]` - post_df_liz$`prob[20]`
quantile(diff, 0.025)
quantile(diff, 0.975)

##### Section 3 Question 10 #####