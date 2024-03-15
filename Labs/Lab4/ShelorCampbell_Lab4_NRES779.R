##Problem 1)
##Load in the necessary packages
library(plyr)
##Set seed to 2021 for standardized results
set.seed(2021)
theta = 6.4
##Put in the formula for simulating a Poisson distribution
y=rpois(50,theta)
y ##Just to double check things
var(y)


##Plot the damn thing
hist(y)
##Use Plyr package to get more fine grained data
plot(count(y), type='h', col=4, lwd =4)

##Set prior values
mu.prior<-10.2
sigma.prior<-0.5
step<-0.01
theta=seq(0,15,step)


##Gamma Distribution - Get Prior Distribution
prior=function(theta, mu=10.2, sigma=0.5){
  shape=mu^2/sigma^2
  rate=mu/sigma^2
  dgamma(theta, shape=shape, rate=rate, log=FALSE)
}
x=theta
yprior=prior(theta,mu.prior,sigma.prior)
plot(theta,yprior, type='l')

##Gamma Distribution - 100000 Simulations
G1<-rgamma(100000, shape=mu.prior^2/sigma.prior^2,
           rate=mu.prior/sigma.prior^2)
mean(G1)
sd(G1)


##Coding a Likelihood Function
like = function(theta,y){
  L=rep(0,length(theta))
  for(i in 1:length(theta))
    L[i]=prod(dpois(y, theta[i],log=FALSE))
  return(L)
}

##Plot the Likelihood Function
plot(theta, like(theta,y=y), type="l", xlim=c(5,15),
     main="Likelihood", xlab = expression(theta), 
     ylab = expression(paste("[y|",theta,"]")))



##Create a Function of the Joint Distribution
joint=function(theta){
  like(theta, y=y)*prior(theta,mu.prior,sigma.prior)
}

##Plot the Joint Distribution
plot(theta, joint(theta), type="l", xlim=c(5,15),
     main="Joint", xlab = expression(theta), 
     ylab = expression(paste("[y|",theta,"] x [",theta,"]")))



##Get Marginal Probability
Marg_Probability<-sum(like(theta,y)*prior(theta,mu.prior,sigma.prior)*step)
Marg_Probability



##Divide Joint Distribution over the Integral of the Joint Distribution
posterior=joint(theta)/Marg_Probability

##Plot the Posterior as a Function of Theta
plot(theta, posterior, type='l', xlim=c(5,15),
     main="Posterior", xlab=expression(theta),
     ylab=expression(paste("[",theta,"|y}")))


##Create a 2x3 Plot of prior graphs
par(mfrow=c(2,3))
##Prior Graph
plot(theta, prior(theta), type='l', main="Prior", xlim=c(5,15),
     xlab=expression(theta),
     ylab=expression(paste("[",theta,"]")))
##Histogram
hist(y, freq=FALSE, breaks=10,main='Histogram of Data')
##Histogram but Skinny
plot(count(y), type="h")
##likelihood Plot
plot(theta, like(theta,y=y), type='l', main="Likelihood", xlim=c(5,15),
     xlab=expression(theta),
     ylab=expression(paste("[y|",theta,"]")))
##Joint Distribution Plot
plot(theta, joint(theta), type='l', main="Joint", xlim=c(5,15),
     xlab=expression(theta),
     ylab=expression(paste("[y|",theta,"] x [",theta,"|y]")))
##Posterior Graph
plot(theta, posterior, type='l', main="Posterior",
     xlim=c(5,15),
     xlab=expression(theta),
     ylab=expression(paste("[",theta,"}y]")))
##Reset Graph Size
par(mfrow=c(1,1))


##Divided Max Posterior by Max Likelihood to Re-Scale Likelihood
(c<-max(posterior)/max(like(theta,y)))
rescaled_like=c*like(theta,y)


##Overlay the plots on one another
plot(theta, rescaled_like,type="l", col='red', main='Scaled Overlay', xlim=c(5,15),
     xlab=expression(theta),ylab="Probability Density")
lines(theta, posterior,type='l')
lines(theta,prior(theta),col='blue')
legend(12,1,c("Scaled Likelihood","Posterior","Prior"),
       lty=rep('solid',3),
       col=c('red','black','blue'))


##Calculate Posterior using Poisson-Gamma Distribution
##Overlay Results
conjugate<-dgamma(theta, sum(y)+mu.prior^2/sigma.prior^2,
                  length(y)+mu.prior/sigma.prior^2)
par(mfrow=c(1,1))
plot(theta,rescaled_like,type='l',main="Scaled Overlay",col='red',xlim=c(5,15),
     xlab=expression(theta),ylab="Probability Density")
lines(theta, prior(theta),col='blue')
lines(theta,conjugate,type='l',lwd=1,col='orange')
lines(theta,posterior,type='l',lwd=1,col='black',lty='dashed')
legend(11,1.2,legend=c("Scaled Likelihood","Prior","Integrated Posterior",
                       "Conjugate Posterior"),
       cex=0.8, lwd=2, bty='n', col=c("red",'blue','orange','black'),
       lty=c("solid",'solid','dashed','solid'))