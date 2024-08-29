##### NRES 779 Lab 6 #####
### Griffin Shelor
### much code (starting from the "Example" section and below) written by Hannah Potts

library(pacman)
p_load(tidyverse, here, invgamma)

##### The Metropolis-Hastings Step #####
## Question 1

set.seed(1874)

# Defining Values
mu_z <- 12
sigma_z <- 1.5
z <- 10

mu_xi <- 10
sigma_xi <- 1.5
xi <- 12

# Normal Distribution Probability Densities
pdf_z <- dnorm(z, mean = mu_z, sd = sigma_z)
pdf_xi <- dnorm(xi, mean = mu_xi, sd = sigma_xi)

print(pdf_z) #0.10934
print(pdf_xi) #0.10934

# Moment Matching for Gamma
variance_1 <- sigma_z^2
shape_1 <- mu_z^2 / variance_1
rate_1 <- mu_z / variance_1

variance_2 <- sigma_xi^2
shape_2 <- mu_xi^2 / variance_2
rate_2 <- mu_xi / variance_2

# Gamma Distribution Probability Densities
gamma_pdf_z <- dgamma(z, shape = shape_1, rate = rate_1)
gamma_pdf_xi <- dgamma(xi, shape = shape_2, rate = rate_2)

print(gamma_pdf_z) #0.1170333
print(gamma_pdf_xi) #0.1008311

##### Deterministic Model #####
g <- function(alpha, gamma, c, L){
  mu <- alpha * (L - c) / (alpha/gamma + (L - c))
  return(mu)
}

##### Data Simulation #####
getData <- function(alpha, gamma, c, sigma){
  set.seed(4)
  par(mfrow=c(1, 1))
  L = sort(runif(50, min = 12, max = 100))
  mu = g(alpha = alpha, gamma = gamma, c = c, L = L)
  y = rgamma(length(mu), mu^2/sigma^2, mu/sigma^2)
  plot(L, y, pch = 16, cex = 1.4, col = "blue")
  lines(L, mu, lwd = 2)
  
  model = nls(y ~ g( alpha = alpha, gamma = gamma, c = c, L = L),
              start = list(alpha = 50, gamma = 4, c = 2))
  
  s = summary(model)
  p = coef(model)
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

setup <- function(n.iter, n.chain, parameter.names, dim.x){
  # Set up storage for chains
  x = list()
  for(i in 1:length(parameter.names)){
    x[[i]] = array(NA, dim = c(dim.x[i], n.iter, n.chain))
  }
  # Assign parameter names to elements of the list
  names(x) = parameter.names
  
  # Enter initial values for parameters here - based on what we had done on Lab 3, I chose 50 for alpha, 2 for gamma, 3 for c, and 12 for sigma 
  x$alpha[1, 1, 1] = 60
  x$c[1, 1, 1] = 10
  x$gamma[1, 1, 1] = 3
  x$sigma[1, 1, 1] = 5
  
  # Enter tuning parameters here - again based on what we had done for lab 3, I chose 5 for alpha, 0.5 for gamma, 1 for c, and 2 for sigma 
  tune=list(
    alpha=10,
    c = 10,
    gamma=0.3,
    sigma=5
  )
  x$tune = tune
  return(x)
}

#just set this to 10 iterations

x <- setup(n.iter = 10, n.chain = 1,
          parameter.names = c("alpha", "c", "gamma", "sigma",  "growth_ratio"), dim.x = c(1, 1, 1, 1, 1))


##### Proposal Function #####
q <- function(theta, mu, tune, type){
  sigma = tune
  if (type == "density") return (dgamma(theta, mu^2/sigma^2, mu/sigma^2))
  if (type == "draw") return (rgamma(1, mu^2/sigma^2, mu/sigma^2))
}

##### Prior Function #####
prior <- function(param, theta){
  if(param == "alpha") return(dunif(theta, min = 0, max = 500, log = TRUE))
  if(param == "c") return(dunif(theta, min = 0, max = 200, log = TRUE))
  if(param == "gamma") return(dunif(theta, min = 0, max = 200, log = TRUE))
  if(param == "sigma" ) return(dunif(theta, min = 0, max = 200, log = TRUE))
}


##### Likelihood Function #####
Like <- function(y, L, alpha, gamma, c, sigma){
  mu = g(alpha = alpha, gamma = gamma, c = c, L = L)
  LogL = dnorm(y, mean = mu, sd = sigma, log = TRUE)
  return(sum(LogL[is.finite(LogL)]))
}

##### Choose Function #####
choose <- function(x, z, Like_z, Like_x, param, tune){
  # These are both logs so we add rather than multiply
  numerator = Like_z + prior(param, theta = z)
  denominator = Like_x + prior(param, theta = x)
  q.ratio = q(theta = x, mu = z, tune = tune, type = "density") / q(theta = z, mu = x, tune = tune, type = "density")
  ### Because these are logs, we exponentiate the difference to get the ratio.
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
  q.y.hat = apply(x$y.hat[,burnin:n.iter, 1], 1,
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

##### Example #####
#generating data
data = getData(alpha = 38.5, gamma = 1.7, sigma = 2, c = 8)

#setting up
n.iter = 50000
x = setup(n.iter = n.iter, n.chain = 1,
          parameter.names = c("alpha", "c", "gamma", "sigma",
                              "y.hat", "growth_ratio"), dim.x = c(1, 1, 1, 1, 50, 1))

#defining parameters
burnin = 15000
tune = x$tune


for (i in 2:50000) {
  #alpha
  a_proposal = q(50, x$alpha[i-1], tune$alpha, "draw")
  
  Like_z = Like(data$y, data$x, alpha=a_proposal, c=x$c[i-1], gamma=x$gamma[i-1], sigma=x$sigma[i-1])
  Like_x = Like(data$y, data$x, alpha=x$alpha[i-1], c=x$c[i-1], gamma=x$gamma[i-1], sigma=x$sigma[i-1])
  chosen = choose(x$alpha[i-1], a_proposal, Like_z, Like_x, "alpha", tune$alpha)
  
  x$alpha[i] = chosen
  
  #c
  c_proposal = q(3, x$c[i-1], tune$c, "draw")
  
  Like_z = Like(data$y, data$x, alpha=x$alpha[i], c=c_proposal, gamma=x$gamma[i-1], sigma=x$sigma[i-1])
  Like_x = Like(data$y, data$x, alpha=x$alpha[i], c=x$c[i-1], gamma=x$gamma[i-1], sigma=x$sigma[i-1])
  chosen = choose(x$c[i-1], c_proposal, Like_z, Like_x, "c", tune$c)
  
  x$c[i] = chosen
  
  #gamma
  g_proposal = q(2, x$gamma[i-1], tune$gamma, "draw")
  
  Like_z = Like(data$y, data$x, alpha=x$alpha[i], c=x$c[i], gamma=g_proposal, sigma=x$sigma[i-1])
  Like_x = Like(data$y, data$x, alpha=x$alpha[i], c=x$c[i], gamma=x$gamma[i-1], sigma=x$sigma[i-1])
  chosen = choose(x$gamma[i-1], g_proposal, Like_z, Like_x, "gamma", tune$gamma)
  
  x$gamma[i] = chosen
  
  #sigma
  s_proposal = q(12, x$sigma[i-1], tune$sigma, "draw")
  
  Like_z = Like(data$y, data$x, alpha=x$alpha[i], c=x$c[i], gamma=x$gamma[i], sigma=s_proposal)
  Like_x = Like(data$y, data$x, alpha=x$alpha[i], c=x$c[i], gamma=x$gamma[i], sigma=x$sigma[i-1])
  chosen = choose(x$sigma[i-1], s_proposal, Like_z, Like_x, "sigma", tune$sigma)
  
  x$sigma[i] = chosen
  
  #y.hat
  
  x$y.hat[,i,1] = g(alpha= x$alpha[i], gamma = x$gamma[i], c= x$c[i], L=data$x)
  
  #growth ratio
  
  x$growth_ratio[i] = x$alpha[i]/x$gamma[i]
}

### plotting
tracePlots(data = data, x = x, n.iter = n.iter, burnin = burnin)
predictPlots(data = data, x = x, n.iter = n.iter, burnin = burnin)
densityPlots(data = data, x = x, n.iter = n.iter, burnin = burnin)

##### Tuning #####
class(x)

tuner = x$alpha[,,]
length(unique(tuner))
length(unique(tuner))/50000

##### run MCMC #####
run_MCMC <- function(alpha, gamma, sigma, c, n.iter,  burnin, atune) {
  
  n.iter = n.iter
  burnin = burnin
  
  # Generating data
  data <- getData(alpha = alpha, gamma = gamma, sigma = sigma, c = c)
  
  # Setting up
  x <- setup(n.iter = n.iter, n.chain = 1,
             parameter.names = c("alpha", "c", "gamma", "sigma", "y.hat", "growth_ratio"), dim.x = c(1, 1, 1, 1, 50, 1))
  
  # Define tuning parameters
  tune <- x$tune
  tune$alpha <- atune
  
  for (i in 2:n.iter) {
    # Alpha
    a_proposal <- q(50, x$alpha[i-1], tune$alpha, "draw")
    
    Like_z <- Like(data$y, data$x, alpha = a_proposal, c = x$c[i-1], gamma = x$gamma[i-1], sigma = x$sigma[i-1])
    Like_x <- Like(data$y, data$x, alpha = x$alpha[i-1], c = x$c[i-1], gamma = x$gamma[i-1], sigma = x$sigma[i-1])
    chosen <- choose(x$alpha[i-1], a_proposal, Like_z, Like_x, "alpha", tune$alpha)
    
    x$alpha[i] <- chosen
    
    # c
    c_proposal <- q(3, x$c[i-1], tune$c, "draw")
    
    Like_z <- Like(data$y, data$x, alpha = x$alpha[i], c = c_proposal, gamma = x$gamma[i-1], sigma = x$sigma[i-1])
    Like_x <- Like(data$y, data$x, alpha = x$alpha[i], c = x$c[i-1], gamma = x$gamma[i-1], sigma = x$sigma[i-1])
    chosen <- choose(x$c[i-1], c_proposal, Like_z, Like_x, "c", tune$c)
    
    x$c[i] <- chosen
    
    # gamma
    g_proposal <- q(2, x$gamma[i-1], tune$gamma, "draw")
    
    Like_z <- Like(data$y, data$x, alpha = x$alpha[i], c = x$c[i], gamma = g_proposal, sigma = x$sigma[i-1])
    Like_x <- Like(data$y, data$x, alpha = x$alpha[i], c = x$c[i], gamma = x$gamma[i-1], sigma = x$sigma[i-1])
    chosen <- choose(x$gamma[i-1], g_proposal, Like_z, Like_x, "gamma", tune$gamma)
    
    x$gamma[i] <- chosen
    
    # sigma
    s_proposal <- q(12, x$sigma[i-1], tune$sigma, "draw")
    
    Like_z <- Like(data$y, data$x, alpha = x$alpha[i], c = x$c[i], gamma = x$gamma[i], sigma = s_proposal)
    Like_x <- Like(data$y, data$x, alpha = x$alpha[i], c = x$c[i], gamma = x$gamma[i], sigma = x$sigma[i-1])
    chosen <- choose(x$sigma[i-1], s_proposal, Like_z, Like_x, "sigma", tune$sigma)
    
    x$sigma[i] <- chosen
    
    # y.hat
    x$y.hat[, i, 1] <- g(alpha = x$alpha[i], gamma = x$gamma[i], c = x$c[i], L = data$x)
    
    # growth ratio
    x$growth_ratio[i] <- x$alpha[i] / x$gamma[i]
  }
  
  tracePlots = function(data, x, n.iter, burnin){
    plot(x$alpha[1,burnin:n.iter,1], type = "l", xlab = "Iteration",
         col = "chartreuse4", ylab = expression(alpha), xlim =c(0,100))
    abline(h = mean(x$alpha[1,burnin:n.iter,1]), col = "red")
  }
  
  # Plotting
  trace <- tracePlots(data = data, x = x, n.iter = n.iter, burnin = burnin)
  #predict <- predictPlots(data = data, x = x, n.iter = n.iter, burnin = burnin)
  #density <- densityPlots(data = data, x = x, n.iter = n.iter, burnin = burnin)
  
  # Return both MCMC results and plots
  return(list(x = x, plots = list(trace = trace)))
  #return(list(x = x, plots = list(trace = trace, predict = predict, density = density))
}

# Example usage
result <- run_MCMC(alpha = 38.5, gamma = 1.7, sigma = 2, c = 8, n.iter = 50000, burnin = 15000, atune = 10)


##### Testing various tuning parameters #####
#2
result <- run_MCMC(alpha = 38.5, gamma = 1.7, sigma = 2, c = 8, n.iter = 50000, burnin = 15000, atune = 2)

#10 - original
result <- run_MCMC(alpha = 38.5, gamma = 1.7, sigma = 2, c = 8, n.iter = 50000, burnin = 15000, atune = 10)

#20
result <- run_MCMC(alpha = 38.5, gamma = 1.7, sigma = 2, c = 8, n.iter = 50000, burnin = 15000, atune = 20)

#50
result <- run_MCMC(alpha = 38.5, gamma = 1.7, sigma = 2, c = 8, n.iter = 50000, burnin = 15000, atune = 50)