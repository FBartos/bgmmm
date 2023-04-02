n1     <- 750
n2     <- 750

mu1    <- rep( -.50, 10)
mu2    <- rep(  .50, 10)

sigma1 <- matrix( 0,  ncol = length(mu1), nrow = length(mu1))
sigma2 <- matrix( .3, ncol = length(mu2), nrow = length(mu2))

diag(sigma1) <- 1
diag(sigma2) <- 1

y1 <- MASS::mvrnorm(n1, mu1, sigma1)
y2 <- MASS::mvrnorm(n2, mu2, sigma2)
y <- rbind(y1, y2)


fit_1500 <- bgmmm(y, 8, parallel = F)

fit_2000
fit_1500
fit_1000
fit_500

fit_2000$models$MC_2
fit_1500$models$MC_2
fit_1000$models$MC_2
fit_500$models$MC_2

EPDS <- readRDS("C:/Projects/papers/networks mixture/t3.EPDS.RDS")
EPDS <- na.omit(EPDS)
apply(EPDS, 2, table)
fit_EPDS <- bgmmm(EPDS, 5, parallel = T)
fit_EPDS$models$MC_2
fit_EPDS$models$MC_4
fit_EPDS$inference$models_marglik

library(LaplacesDemon)
R <- 10
data <- list(
  prior_k = R + 1,
  prior_R = diag(1, R),
  R       = R
)

cors <- sapply(1:1000, function(i){
  Sigma <- Prec2Cov(LaplacesDemon::rwishart(data$prior_k, data$prior_R))
  Cors  <- cov2cor(Sigma)
  Cors[upper.tri(Cors, F)]
}, simplify = F)
hist(unlist(cors), breaks = 100, xlim = c(-1, 1))

sds <- unlist(sapply(1:10000, function(i){
  Sigma <- Prec2Cov(LaplacesDemon::rwishart(data$prior_k, data$prior_R))
  sqrt(diag(Sigma))
}, simplify = F))
hist(sds[sds<10],   breaks = 100)
hist(sds^2, breaks = 100)

min(sds)
min(sds^2)
mean(sds)
mean(sds^2)
fit1$models$MC_2

fit1$inference
round(fit$inference$post_prob$K, 3)
round(fit$inference$post_prob$model_types, 3)
fit$models[[8]]$posterior$mu

f_mc  <- fit_jags(y, 1, "mc")
f_MC2 <- fit_jags(y, 2, "MC")
f_mC2 <- fit_jags(y, 2, "mC")
f_Mc2 <- fit_jags(y, 2, "Mc")
f_mc$marg_lik


R <- 10
log(2^R) + sum( sapply(1:R,function(x)log( chol(Cov2Prec(sigma1))[x,x]^(R-x+2)) ) )
log(2^R) + sum( sapply(1:R,function(x)log( chol(Cov2Prec(sigma2))[x,x]^(R-x+2)) ) )



cors <- unlist(sapply(f_mc$posterior$Sigma, function(x)cov2cor(x)[upper.tri(x)]))

hist(unlist(sapply(f_MC2$posterior$Sigma[[1]], function(x)cov2cor(x)[upper.tri(x)])))
hist(unlist(sapply(f_MC2$posterior$Sigma[[2]], function(x)cov2cor(x)[upper.tri(x)])))

hist(unlist(sapply(f_mC2$posterior$Sigma[[1]], function(x)cov2cor(x)[upper.tri(x)])))
hist(unlist(sapply(f_mC2$posterior$Sigma[[2]], function(x)cov2cor(x)[upper.tri(x)])))


round(bridgesampling::post_prob(f_mc$marg_lik, f_Mc2$marg_lik, f_mC2$marg_lik, f_MC2$marg_lik), 100)
round(bridgesampling::post_prob(f_Mc2$marg_lik, f_mC2$marg_lik, f_MC2$marg_lik), 3)

chol(prec2)
chol(sigma2)
chol2inv(chol(sigma2))
prec2 <- Cov2Prec(sigma2)
Prec2Cov(t(chol(prec2)) %*% chol(prec2))
