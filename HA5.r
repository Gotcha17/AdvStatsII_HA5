##### a) ##################################################################################################

# Log-likelihood function for the normal dirstribution
# mu and sigma are inputed through the "par" vector, this is needed due to the use of the optim() fct.
loglike = function(par, x){
  mu = par[1]
  sigma = par[2]
  ll = -length(x)/2*log(2*pi*sigma^2)-1/(2*sigma^2)*sum(x^2)+mu/(sigma^2)*sum(x)-length(x)/(2*sigma^2)*mu^2
  return(ll)
}

##### b) ##################################################################################################

# Compute moment function for N
moments <- function(par, x){
  mu <- par[1]
  sigma <- par[2]
  m1 <- x - mu                       # first noncentral moment
  m2 <- x^2 - mu^2 - sigma^2         # second noncentral moment
  return(cbind(m1, m2))              # have to return a matrix for using it in gmm-function
}

##### c) ##################################################################################################

# Construct sample of size 100 from N(5,1)-distribution
set.seed(666)
sample <- rnorm(100, mean = 5, sd = 2)

##### d) ##################################################################################################

# Compute LS estimator
lm.fit <- lm(sample ~ 1)
lm.fit <- c(lm.fit$coefficients, sigma(lm.fit))
names(lm.fit) <- c("mu", "sigma")
lm.fit

# Compute ML estimator
ml.fit <- optim(par=c(0,1), fn=loglike, x = sample, control=list(fnscale=-1))
ml.fit <- ml.fit$par
names(ml.fit) <- c("mu", "sigma")
ml.fit

# Compute GMM estimator

#install.packages("gmm")
library("gmm")
gmm.fit <- gmm(g = moments, x = sample, t0 = c(0, 1))
gmm.fit <- gmm.fit$coefficients
names(gmm.fit) <- c("mu", "sigma")
gmm.fit

##### e) ##################################################################################################
MSE <- matrix(nrow = 3, ncol = 2, dimnames = list(c("MSE LS", "MSE ML", "MSE GMM"),c("mu", "sigma")))
# MSE LS
MSE[1,1] <- (lm.fit[1]-5)^2
MSE[1,2] <- (lm.fit[2]-2)^2

# MSE ML
MSE[2,1] <- (ml.fit[1]-5)^2
MSE[2,2] <- (ml.fit[2]-2)^2

# MSE GMM
MSE[3,1] <- (gmm.fit[1]-5)^2
MSE[3,2] <- (gmm.fit[2]-2)^2
MSE

# As we can see from the MSE matrix, we achieve lowest MSE for estimator of mu and sigma with the GMM estimation
###########################################################################################################
