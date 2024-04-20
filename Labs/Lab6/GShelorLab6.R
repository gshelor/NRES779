## NRES 779 Lab 6
## Griffin Shelor

library(pacman)
p_load(tidyverse, here, invgamma)

##### The Metropolis-Hastings Step #####
## Question 1

set.seed(1874)
prob_dens_z <- dnorm(10, 12, 1.5)
prob_dens_z

prob_dens_xi <- dnorm(12, 10, 1.5)
prob_dens_xi

## using gamma distribution
# moment matching

##### Deterministic Model #####
g <- function(alpha, gamma, c, L){
  mu <- alpha * (L - c) / (alpha/gamma + (L - c))
  return(mu)
}

##### Data Simulation #####
getData = function(alpha, gamma, c, sigma){
  set.seed(4)
  par(mfrow = c(1, 1))
  L <- sort(runif(50, min = 12, max = 100))
  mu <- g(alpha = alpha, gamma = gamma, c = c, L = L)
  y <- rgamma(length(mu), mu^2/sigma^2, mu/sigma^2)
  plot(L, y, pch = 16, cex = 1.4, col = "blue")
  lines(L, mu, lwd = 2)
  model <- nls(y ~ g( alpha = alpha, gamma = gamma, c = c, L = L),
              start = list(alpha = 50, gamma = 4, c = 2))
  s <- summary(model)
  p <- coef(model)
  y.hat = g(alpha = p[1], gamma = p[2], c = p[3], L = L)
  lines(L, y.hat, col = "red", lwd = 2)
  legend(40, 18, c("Generating values", "NLS fit"),
         lty = c("solid", "solid"), col = c("black", "red"),
         bty = "n", lwd = 2)
  return(list(
    x = L,
    y = y,
    nls.alpha = p[1],
    nls.gamma = p[2],
    nls.c = p[3],
    nls.sigma = s$sigma,
    gen.alpha = alpha,
    gen.gamma = gamma,
    gen.c = c,
    gen.sigma = sigma))
}


##### Setup #####

setup = function(n.iter, n.chain, parameter.names, dim.x){
  # Set up storage for chains
  x <- list()
  for(i in 1:length(parameter.names)){
    x[[i]] <- array(NA, dim = c(dim.x[i], n.iter, n.chain))
  }
  # Assign parameter names to elements of the list
  names(x) <- parameter.names
  # Enter initial values for parameters here
  x$alpha[1, 1, 1] = 60
  x$c[1, 1, 1] = 10
  x$gamma[1, 1, 1] = 3
  x$sigma[1, 1, 1] = 5
  # Enter tuning parameters here
  tune=list(
    alpha=10,
    c = 1,
    gamma=.3,
    sigma=2
  )
  x$tune = tune
  return(x)
}

## calling setup function
# x = setup(n.iter = n.iter, n.chain = 1,
#           parameter.names = c("alpha", "c", "gamma", "sigma",
#                               "y.hat", "growth_ratio"), dim.x = c(1, 1, 1, 1, n, 1))

##### Proposal Function #####

q <- function(theta, mu, tune, type){
  sigma = tune
  if (type == "density") return (dgamma(theta, muˆ2/sigmaˆ2, mu/sigmaˆ2))
  if (type == "draw") return (rgamma(1, muˆ2/sigmaˆ2, mu/sigmaˆ2))
}

##### Prior Function #####

prior <- function(param, theta){
  if(param == "alpha") return( dunif(theta, min = 0, max = 500, log = TRUE))
  if(param == "c") return(dunif(theta, min = 0, max = 200, log = TRUE))
  if(param == "gamma") return(dunif(theta, min = 0, max = 200, log = TRUE))
  if(param == "sigma" ) return(dunif(theta, min = 0, max = 200, log = TRUE))
}

##### Likelihood Function #####

Like <- function(y, L, alpha, gamma,c ,sigma){
  mu=g(alpha = alpha, gamma = gamma, c = c, L = L)
  LogL = dnorm(y, mean = mu, sd = sigma, log = TRUE)
  return(sum(LogL[is.finite(LogL)]))
}

##### Choose Function #####

choose <- function(x, z, Like_z, Like_x, param, tune){
  # These are both logs so we add rather than multiply
  numerator = Like_z + prior(param, theta = z)
  denominator = Like_x + prior(param, theta = x)
  q.ratio = q(theta = x, mu = z, tune = tune, type = "density") /
    q(theta = z, mu = x, tune = tune, type = "density")
  # Because these are logs, we exponentiate the difference to get the ratio.
  R = exp(numerator - denominator) * q.ratio
  if (R > runif(1, min = 0, max = 1)) new = z else new = x
  return(new)
}

## example calling the function
# choose(x = x$gamma[i - 1], z = z, Like_z = Like_z, Like_x = Like_x,
#        param = "gamma", tune = tune$gamma)

##### Plotting Functions #####

tracePlots <- function(data, x, n.iter, burnin){
  par(mfrow = c(2, 2))
  plot(x$gamma[1,burnin:n.iter,1], type = "l", xlab = "Iteration",
       col = "chartreuse4", ylab = expression(gamma))
  abline(h = mean(x$gamma[1,burnin:n.iter,1]), col = "red")
  plot(x$alpha[1,burnin:n.iter,1], type = "l", xlab = "Iteration",
       col = "chartreuse4", ylab = expression(alpha))
  abline(h = mean(x$alpha[1,burnin:n.iter,1]), col = "red")
  plot(x$c[1,burnin:n.iter,1], type = "l", xlab = "Iteration",
       col = "chartreuse4", ylab = expression(c))
  abline(h = mean(x$c[1,burnin:n.iter,1]), col = "red")
  plot(x$sigma[1,burnin:n.iter,1], type = "l", xlab = "Iteration",
       col = "chartreuse4", ylab = expression(sigma))
  abline(h = mean(x$sigma[1,burnin:n.iter,1]), col="red")
}
predictPlots <- function(data, x, n.iter, burnin){
  par(mfrow = c(1, 1))
  q.y.hat = apply(x$y.hat[, burnin:n.iter, 1], 1,
                  function(x) quantile(x, c(.025, .5, .975)))
  plot(data$x, data$y, xlab = "Light level",
       ylab = "Growth rate", main = "Prediction of growth rate",
       pch = 16, cex = 1.4, col = "blue")
  lines(data$x, g(alpha = data$nls.alpha,
                  gamma = data$nls.gamma, c = data$nls.c, L = data$x),
        col="blue", lwd = 4)
  lines(data$x,q.y.hat[2,], col = "orange", lwd = 2)
  lines(data$x,q.y.hat[1,], lty = "dashed", col = "orange", lwd = 2)
  lines(data$x,q.y.hat[3,], lty = "dashed", col = "orange", lwd = 2)
  legend(40,18, c("Median", "2.5% quantile",
                  "97.5% quantile", "NLS fit"),
         lty = c("solid", "dashed", "dashed"),
         col = c("orange", "orange", "orange", "blue"),
         lwd = c(2, 2, 2, 4),
         bty = "n")
}
plot_density <- function(p, v1, v2, param, burnin, n.iter){
  hist(p[burnin:n.iter], breaks = 40, xlab = "Value of parameter",
       freq = FALSE, main = param, col = "gray")
  abline(v = v1, col = "red", lwd = 3)
  abline(v = v2, col = "blue", lwd = 2)
  abline(v = median(p[burnin:n.iter]), col = "orange", lwd = 2)
}
densityPlots <- function(data, x, n.iter, burnin){
  par(mfrow = c(2, 2))
  plot_density(p = x$alpha, v1 = data$gen.alpha,
               v2 = data$nls.alpha, param = expression(alpha),
               burnin = burnin, n.iter = n.iter)
  plot_density(p = x$gamma, v1 = data$gen.gamma,
               v2 = data$nls.gamma, param = expression(gamma),
               burnin = burnin, n.iter = n.iter)
  plot_density(p = x$c, v1 = data$gen.c, v2 = data$nls.c,
               param = expression(c), burnin = burnin,
               n.iter = n.iter)
  plot_density(p = x$sigma, v1 = data$gen.sigma, v2 = data$nls.sigma,
               param = expression(sigma),
               burnin = burnin, n.iter = n.iter)
}

##### Problem #####
## setting number of iterations
n_iter = 50000

data = getData(50, 3, 2, 12)

