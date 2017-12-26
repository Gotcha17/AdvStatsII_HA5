
# Problem 1 ---------------------------------------------------------------

## a)

# define the log-likelihood function
logilog3 <- function(param, x){
 ll <- -length(x)/2*log(2*pi*param[2]^2)-1/(2*param[2]^2)*sum(x^2)+param[1]/param[2]^2*sum(x)-length(x)/2*(param[1]/param[2])^2
 return(ll)
}

## b)

#define a function to compute the euclidean length of a vector
#plug in difference of two vectors to obtain their euclidean distance
disti <- function(x){
 d <- sqrt(sum(x^2))
 return(d)
}

## c) - d)

# define global variables
n <- seq(from = 100, to = 1000, by = 100)                  # sample sizes
mu_theo <- seq(from = -5, to = 5, length.out = 100)        # values for mu
sigma_theo <- seq(from = 0.1, to = 25, length.out = 100)   # values for sigma
#create array to store distance for all possible combinations of n, mu and sigma
dist <- array(NA, c(length(n), length(mu_theo), length(sigma_theo)))

start <- proc.time()
row <- 1
for (i in n){
 col <- 1
 for (j in mu_theo){
   col2 <- 1
   for (k in sigma_theo){
     x <- rnorm(i, mean = j, sd = k)
     #calculate ML estimates for mu and sigma (with starting values 0 and 1)
     emp_par <- optim(par = c(0,1), fn = logilog3, x = x , control = list(fnscale = -1))$par
     theo_par <- c(j, k)
     #calculate euclidean distance between empirical and true parameters
     dist[row, col, col2] <- disti(emp_par-theo_par)
     col2 <- col2 + 1
   }
 col <- col + 1
}
row <- row + 1
}
proc.time()-start

## e)

#use different starting values for mu and sigma to show that the results
#from (d) depend on these values
sample <- rnorm(100, mean = 4, sd = 2)
#starting values for mu:
mu_start <- seq(from = -1000, to = 1000, length.out = 10)
#starting values for sigma:
sigma_start <- seq(from = 0.01, to = 250, length.out = 10)
#matrix to store parameters for different starting value combinations:
 data <- matrix(NA, nrow = length(mu_start)*length(sigma_start), ncol = 2)
row <- 1
for (i in mu_start){
 for (j in sigma_start){
   emp_par <- optim(par = c(i,j), fn = logilog3, x = sample , control = list(fnscale = -1))$par
 data[row, 1] <- emp_par[1]
   data[row, 2] <- emp_par[2]
   row <- row + 1
 }
}

print(data[1:40,])

## For sigma=0.01 as starting value for sigma the estimates are really bad (see lines 1,11,21,...,91
# of the data matrix). For all other starting value combinations the estimates are quite simi-
# lar, only the sign of the estimate for sigma alternates (which means that the variance is esti-
## mated correctly).


# Problem 2 ---------------------------------------------------------------

## a)

# compute the log-likelihood function of N(mu,1)-distributed sample
log.norm <- function(mu, x){
 ll <- length(x)*log(1/sqrt(2*pi))-1/2*sum(x^2)+mu*sum(x)-length(x)/2*mu^2
 return(ll)
}

## b)

# Compute moment function
moments <- function(mu, x){
 m1 <- x - mu                 # first noncentral moment
 m2 <- x^2 - mu^2 - 1         # second noncentral moment
 return(cbind(m1, m2))       # have to return a matrix for using it in gmm-function
}

## c)

# Construct sample of size 100 from N(5,1)-distribution
set.seed(666)
sample <- rnorm(100, mean = 5, sd = 1)

## d)

# Compute LS estimator
lm.fit <- lm(sample ~ 1)
lm.fit$coefficients

# Compute ML estimator
ml.fit <- optimize(f = log.norm, lower = -10000, upper = 10000, x = sample, maximum = T)
ml.fit$maximum

# Compute GMM estimator

#install.packages("gmm")
library("gmm")
gmm.fit <- gmm(g = moments, x = sample, t0 = 0, optfct = "optimize", lower = -10000, upper = 10000)
gmm.fit$coefficients

## e)

# All 3 methods used in (d) to calculate an estimate for μ yield the same result.


# Problem 3 ---------------------------------------------------------------

## a)

# Compute log-likelihood function
log.gamma = function(beta, x){
 ll = -sum(x)/beta - length(x)*log(beta)
 return(ll)
}

## b)

# Compute moment function
moments = function(beta, x){
 m1 = x-beta                   #first non-central moment
 m2 = x^2-2*beta^2             #second non-central moment
 return(cbind(m1, m2))
}

## c)

# Construct sample of size 100 from Gamma(1,5)-distribution
sample = rgamma(100, shape = 1, scale = 5)

## d)

# Calculate LS estimator
lm.fit = lm(sample ~ 1)
lm.fit$coefficients

# Calculate ML estimator
ml.fit = optimize(f = log.gamma, lower = 0, upper = 10000, x = sample, maximum = T)
ml.fit$maximum

# Construct GMM estimator
gmm.fit = gmm(g = moments, x = sample, t0 = 0, optfct = "optimize", lower = 0, upper = 10000)
gmm.fit$coefficients

## e)

# Least squares and maximum likelihood yield the same estimate for β. The GMM estimate is slightly larger.