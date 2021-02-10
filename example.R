### setup
# required packages
packages <- c("rjags", "runjags", "MASS", "LaplacesDemon", "coda", "parallel", "devtools")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
# and a special branch of bridgesampling since they didn't fix the bug yet
devtools::install_github("https://github.com/FBartos/bridgesampling/tree/simplex-fix")

### the actual fun
# get some fake data (from my anecdotal evidence ~ 2,000 observations seems to produce quite good outcomes)
# the high amount of parameters needed to estimate the models with differences in covariances probably compounds the issues with priors 
set.seed(666)
K      <- 10
n1     <- 400
n2     <- 600

mu1    <- rep( -.50, K)
mu2    <- rep(  .50, K)

sigma1 <- matrix( 0,  ncol = K, nrow = K)
sigma2 <- matrix( .3, ncol = K, nrow = K)

diag(sigma1) <-  .8
diag(sigma2) <- 1.2

y1 <- MASS::mvrnorm(n1, mu1, sigma1)
y2 <- MASS::mvrnorm(n2, mu2, sigma2)
y  <- rbind(y1, y2)

# load the functions
source('bgmmm_jags.R')

### fit the mixture up to 4 components (with default settings - more in the rest of the code)
# the priors assume that the data are somewhat standardize - with grand mean ~ 0 and grand sd ~ 1 
fit_fake <- bgmmm(y, 4, parallel = TRUE)

# check that the marginals likelihoods were computed (-Inf stands for a failure)
# - using more than 1 chain usually leads to problems in marglik estimation, it might be due to 
#   - different chains have different group labels
#   - the chains do not converge into stationary distribution
#   this also needs to be further investigated, but I did not have time for it yet
fit_fake$inference$models_marglik

# check the main summary
fit_fake

# check the correct model estimates
fit_fake$models$MC_2

### some additional options
# you can change the:
# - number of chains
# - number of iterations for adapting the chains
# - number of iterations for burnin
# - number of sampling iterations
# - maximum number of bridge sampling iterations
# by setting the options argument (using the default values)
options <- list(
  "chains"      = 1,
  "adapt"       = 500,
  "burnin"      = 2000,
  "iter"        = 5000,
  "bridge_iter" = 100000
)
fit_fake <- bgmmm(y, 4, parallel = TRUE, options = options)


### more details regarding the prior distributions
library(LaplacesDemon)

### the prior distribuition for the group proportions is a Dirichlet distribution
# all prior proportions = prior_alpha, defaults to 5% of the data set (as if we already observed 5% observations being in each group)
# - very nice for regularization
# the prior distribution for proportions of participants belonging to each group (with K groups) might be visualized:
K           <- 10
prior_alpha <- nrow(y) * 0.05

hist(rdirichlet(10000, rep(nrow(y) * 0.05, K))[,1], xlim = c(0, 1))

### the prior distribution for means is just multivariate normal with (default is multivariate standard normal)
# all means      = prior_mu, defaults to 0
# all precisions = prior_S,  defaults to 1
# - needs to be re-parametrized for to grand mean and differences from the grand mean for each group
#    - it will also allow many additional options that might be tested (I would explore those only later one), e.g., 
#      - are all differences the same?
#      - are all differences in the same direction?

# the prior distribution for covariance matrices is an inverse-Wishard distribution (equals to a Wishard distribution for precision matricies)
# all degrees of freedom = prior_k, defaults to 2 times the number of columns
# all precisions         = prior_r, defaults to 1 over the number of columns
# - needs to be re-parametrized into two priors:
#   - one for the grand variance matrix
#   - one for the grand correlation matrix
#   and differences from the grand variance/correlation matrix for each group
# the current prior is specified to provide quite a nice shrinkage, but it's terrible for testing
prior_k <- 2*ncol(y)
prior_r <- 1/ncol(y)

# first I sample the covariances and then extract the sds and correlations to create an understandable visualization
prior_covariances <- sapply(1:1000, function(i)Prec2Cov(rwishart(prior_k, diag(prior_r, ncol(y)))), simplify = F)

# prior on correlations
hist(unlist(sapply(prior_covariances, function(cov)Cov2Cor(cov)[lower.tri(cov)])), breaks = 100, xlim = c(-1, 1))

# prior on standard deviations
hist(unlist(sapply(prior_covariances, function(cov)sqrt(diag(cov)))), breaks = 50, xlim = c(0, 5))

### you can fit a model with different prior distributions using the priors argument (using the default values):
priors <- list(
  "prior_alpha" = nrow(y) * 0.05,
  "prior_mu"    = 0,
  "prior_S"     = 1,
  "prior_k"     = 2*ncol(y),
  "prior_r"     = 1/ncol(y)
)
fit_fake <- bgmmm(y, 4, parallel = TRUE, priors = priors)
