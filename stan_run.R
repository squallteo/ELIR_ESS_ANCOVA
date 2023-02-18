library(rstan)
library(mvtnorm)
library(tidyverse)

#total sample size
n_total <- 120
#randomization ratio: treatment:control
RRatio <- 2 
#baseline and post-baseline mean
ped_means_pbo <- c(8.0, 8.0) 
ped_means_active <- c(8.0, 8.0)
#covariance matrix, same for treatment and control arms
ped_vcov <- matrix(c(1.0^2, 0.4, 0.4, 1.43^2), nrow=2)

pbo_dat <- rmvnorm((1/(1+RRatio))*n_total, ped_means_pbo, ped_vcov)
active_dat <- rmvnorm((RRatio/(1+RRatio))*n_total, ped_means_active, ped_vcov)
total_dat <- rbind(cbind(pbo_dat, rep(0, (1/(1+RRatio))*n_total)), cbind(active_dat, rep(1, (RRatio/(1+RRatio))*n_total)))

alldt <- as.data.frame(cbind(total_dat, total_dat[,2]-total_dat[,1]))
colnames(alldt) <- c('Base', "Aval", 'TRTP', 'Chg')
#####################################
prior_m <- -0.6
prior_s <- 0.6

stan_data <- list(N=nrow(alldt), base=alldt$Base, trt=alldt$TRTP,
                  y=alldt$Chg, sigma=sd(alldt$Chg),
                  delta_m=prior_m, delta_s=prior_s)
model <- stan_model("ANCOVA.stan")
fit <- sampling(model, stan_data, chains=3, iter=12000, warmup=2000, thin=4, show_messages=F)

delta_post <- as_tibble(rstan::extract(fit)$delta, .name_repair = "minimal")
